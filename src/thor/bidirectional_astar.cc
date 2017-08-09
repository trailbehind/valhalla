#include <map>
#include <algorithm>
#include "thor/bidirectional_astar.h"
#include "baldr/datetime.h"
#include "midgard/logging.h"

using namespace valhalla::baldr;
using namespace valhalla::sif;

namespace {

// Find a threshold to continue the search - should be based on
// the max edge cost in the adjacency set?
int GetThreshold(const TravelMode mode, const int n) {
  return (mode == TravelMode::kDrive) ?
      n + std::min(8500, std::max(100, n / 3)) :
      n + 500;
}

}

namespace valhalla {
namespace thor {

constexpr uint64_t kInitialEdgeLabelCountBD = 1000000;

// Default constructor
BidirectionalAStar::BidirectionalAStar()
    : PathAlgorithm(),
      access_mode_(kAutoAccess),
      mode_(TravelMode::kDrive),
      travel_type_(0),
      threshold_(0),
      cost_diff_(0.0f) {
}

// Destructor
BidirectionalAStar::~BidirectionalAStar() {
  Clear();
}

// Clear the temporary information generated during path construction.
void BidirectionalAStar::Clear() {
  forward_labels_.clear();
  reverse_labels_.clear();
}

// Initialize the A* heuristic and adjacency lists for both the forward
// and reverse search.
void BidirectionalAStar::Init(const PointLL& origll, const PointLL& destll) {
  // Initialize the A* heuristics
  float factor = costing_->AStarCostFactor();
  astarheuristic_forward_.Init(destll, factor);
  astarheuristic_reverse_.Init(origll, factor);

  // TODO - appropriate sizing for label set elements based on estimated route
  // length. Perhaps some way to recover memory (to work better in low memory
  // environments).

  // Initialize forward LabelSet
  float mincostf  = astarheuristic_forward_.Get(origll);
  forward_labels_.Init(kInitialEdgeLabelCountBD, costing_->UnitSize(),
            kBucketCount, mincostf);

  // Initialize reverse LabelSet
  float mincostr = astarheuristic_reverse_.Get(destll);
  reverse_labels_.Init(kInitialEdgeLabelCountBD, costing_->UnitSize(),
            kBucketCount, mincostr);

  // Set the cost diff between forward and reverse searches (due to distance
  // approximator differences). This is used to "even" the forward and reverse
  // searches.
  cost_diff_ = mincostf - mincostr;

  // Initialize best connection with max cost
  best_connection_ = { GraphId(), GraphId(),
                       std::numeric_limits<float>::max() };

  // Set the threshold to 0 (used to extend search once an initial
  // connection has been found).
  threshold_ = 0;

  // Support for hierarchy transitions
  hierarchy_limits_forward_ = costing_->GetHierarchyLimits();
  hierarchy_limits_reverse_ = costing_->GetHierarchyLimits();
}

// Expand from a node in the forward direction
void BidirectionalAStar::ExpandForward(GraphReader& graphreader,
       const GraphId& node, const BDEdgeLabel& pred, const uint32_t pred_idx,
       const bool from_transition) {
  // Get the tile and the node info. Skip if tile is null (can happen
  // with regional data sets) or if no access at the node.
  const GraphTile* tile = graphreader.GetGraphTile(node);
  if (tile == nullptr) {
    return;
  }
  const NodeInfo* nodeinfo = tile->node(node);
  if (!costing_->Allowed(nodeinfo)) {
    return;
  }

  // Expand from end node in forward direction.
  uint32_t shortcuts = 0;
  GraphId edgeid = { node.tileid(), node.level(), nodeinfo->edge_index() };
  const DirectedEdge* directededge = tile->directededge(edgeid);
  for (uint32_t i = 0; i < nodeinfo->edge_count(); ++i, ++directededge, ++edgeid) {
    // Handle transition edges - expand from the end node of the transition
    // (unless this is called from a transition).
    if (directededge->trans_up()) {
      if (!from_transition) {
        hierarchy_limits_forward_[node.level()].up_transition_count++;
        ExpandForward(graphreader, directededge->endnode(), pred, pred_idx, true);
      }
      continue;
    }
    if (directededge->trans_down()) {
      if (!from_transition &&
          !hierarchy_limits_forward_[directededge->endnode().level()].StopExpanding()) {
        ExpandForward(graphreader, directededge->endnode(), pred, pred_idx, true);
      }
      continue;
    }

    // Quick check to skip if no access for this mode or if edge is
    // superseded by a shortcut edge that was taken.
    if (!(directededge->forwardaccess() & access_mode_) ||
         (shortcuts & directededge->superseded())) {
      continue;
    }

    // Get the current set. Skip this edge if permanently labeled (best
    // path already found to this directed edge).
    EdgeStatusInfo edgestatus = forward_labels_.edgestatus(edgeid);
    if (edgestatus.set() == EdgeSet::kPermanent) {
      continue;
    }

    // Skip this edge if no access is allowed (based on costing method)
    // or if a complex restriction prevents transition onto this edge.
    if (!costing_->Allowed(directededge, pred, tile, edgeid) ||
         costing_->Restricted(directededge, pred, forward_labels_.edgelabels(),
                              tile, edgeid, true)) {
      continue;
    }

    // Get cost. Separate out transition cost. Update the_shortcuts mask.
    // to supersede any regular edge, but only do this once we have stopped
    // expanding on the next lower level (so we can still transition down to
    // that level).
    if (directededge->is_shortcut() &&
        hierarchy_limits_forward_[edgeid.level()+1].StopExpanding()) {
      shortcuts |= directededge->shortcut();
    }
    Cost tc = costing_->TransitionCost(directededge, nodeinfo, pred);
    Cost newcost = pred.cost() + tc + costing_->EdgeCost(directededge);

    // Check if edge is temporarily labeled and this path has less cost. If
    // less cost the predecessor is updated and the sort cost is decremented
    // by the difference in real cost (A* heuristic doesn't change)
    if (edgestatus.set() == EdgeSet::kTemporary) {
      forward_labels_.decrease(edgestatus.index(), pred_idx, newcost, tc);
      continue;
    }

    // Get end node tile (skip if tile is not found) and opposing edge Id
    const GraphTile* t2 = directededge->leaves_tile() ?
        graphreader.GetGraphTile(directededge->endnode()) : tile;
    if (t2 == nullptr) {
      continue;
    }
    GraphId oppedge = t2->GetOpposingEdgeId(directededge);

    // Find the sort cost (with A* heuristic) using the lat,lng at the
    // end node of the directed edge.
    float dist = 0.0f;
    float sortcost = newcost.cost + astarheuristic_forward_.Get(
          t2->node(directededge->endnode())->latlng(), dist);

    // Add edge to the LabelSet
    forward_labels_.add(edgeid, {pred_idx, edgeid, oppedge, directededge,
                  newcost, sortcost, dist, mode_, tc,
                  (pred.not_thru_pruning() || !directededge->not_thru())});
  }
}

// Expand from a node in reverse direction.
void BidirectionalAStar::ExpandReverse(GraphReader& graphreader,
         const GraphId& node, const BDEdgeLabel& pred, const uint32_t pred_idx,
         const DirectedEdge* opp_pred_edge, const bool from_transition) {
  // Get the tile and the node info. Skip if tile is null (can happen
  // with regional data sets) or if no access at the node.
  const GraphTile* tile = graphreader.GetGraphTile(node);
  if (tile == nullptr) {
    return;
  }
  const NodeInfo* nodeinfo = tile->node(node);
  if (!costing_->Allowed(nodeinfo)) {
    return;
  }

  // Expand from end node in reverse direction.
  uint32_t shortcuts = 0;
  GraphId edgeid = { node.tileid(), node.level(), nodeinfo->edge_index() };
  const DirectedEdge* directededge = tile->directededge(edgeid);
  for (uint32_t i = 0; i < nodeinfo->edge_count(); ++i, ++directededge, ++edgeid) {
    // Handle transition edges - expand from the end not of the transition
    // unless this is called from a transition.
    if (directededge->trans_up()) {
      if (!from_transition) {
        hierarchy_limits_reverse_[node.level()].up_transition_count++;
        ExpandReverse(graphreader, directededge->endnode(), pred, pred_idx,
                      opp_pred_edge, true);
      }
      continue;
    } else if (directededge->trans_down()) {
      if (!from_transition &&
          !hierarchy_limits_reverse_[directededge->endnode().level()].StopExpanding()) {
        ExpandReverse(graphreader, directededge->endnode(), pred, pred_idx,
                      opp_pred_edge, true);
      }
      continue;
    }

    // Quick check to skip if no access for this mode or if edge is
    // superseded by a shortcut edge that was taken.
    if (!(directededge->reverseaccess() & access_mode_) ||
         (shortcuts & directededge->superseded())) {
      continue;
    }

    // Get the current set. Skip this edge if permanently labeled (best
    // path already found to this directed edge).
    EdgeStatusInfo edgestatus = reverse_labels_.edgestatus(edgeid);
    if (edgestatus.set() == EdgeSet::kPermanent) {
      continue;
    }

    // Get end node tile, opposing edge Id, and opposing directed edge.
    const GraphTile* t2 = directededge->leaves_tile() ?
        graphreader.GetGraphTile(directededge->endnode()) : tile;
    if (t2 == nullptr) {
      continue;
    }
    GraphId oppedge = t2->GetOpposingEdgeId(directededge);
    const DirectedEdge* opp_edge = t2->directededge(oppedge);

    // Skip this edge if no access is allowed (based on costing method)
    // or if a complex restriction prevents transition onto this edge.
    if (!costing_->AllowedReverse(directededge, pred, opp_edge, t2, oppedge) ||
         costing_->Restricted(directededge, pred, reverse_labels_.edgelabels(),
                              tile, edgeid, false)) {
      continue;
    }

    // Get cost. Use opposing edge for EdgeCost. Update the_shortcuts mask
    // to supersede any regular edge, but only do this once we have stopped
    // expanding on the next lower level (so we can still transition down to
    // that level). Separate the transition seconds so we can properly recover
    // elapsed time on the reverse path.
    if (directededge->is_shortcut() &&
        hierarchy_limits_reverse_[edgeid.level()+1].StopExpanding()) {
      shortcuts |= directededge->shortcut();
    }
    Cost tc = costing_->TransitionCostReverse(directededge->localedgeidx(),
                             nodeinfo, opp_edge, opp_pred_edge);
    Cost newcost = pred.cost() + costing_->EdgeCost(opp_edge);
    newcost.cost += tc.cost;

    // Check if edge is temporarily labeled and this path has less cost. If
    // less cost the predecessor is updated and the sort cost is decremented
    // by the difference in real cost (A* heuristic doesn't change)
    if (edgestatus.set() != EdgeSet::kUnreached) {
      reverse_labels_.decrease(edgestatus.index(), pred_idx, newcost, tc);
      continue;
    }

    // Find the sort cost (with A* heuristic) using the lat,lng at the
    // end node of the directed edge.
    float dist = 0.0f;
    float sortcost = newcost.cost + astarheuristic_reverse_.Get(
       t2->node(directededge->endnode())->latlng(), dist);

    // Add edge to the LabelSet
    reverse_labels_.add(edgeid, {pred_idx, edgeid, oppedge,
        directededge, newcost, sortcost, dist, mode_, tc,
        (pred.not_thru_pruning() || !directededge->not_thru())});
  }
}

// Calculate best path using bi-directional A*. No hierarchies or time
// dependencies are used. Suitable for pedestrian routes (and bicycle?).
std::vector<PathInfo> BidirectionalAStar::GetBestPath(PathLocation& origin,
             PathLocation& destination, GraphReader& graphreader,
             const std::shared_ptr<DynamicCost>* mode_costing,
             const sif::TravelMode mode) {
  // Set the mode and costing
  mode_ = mode;
  costing_ = mode_costing[static_cast<uint32_t>(mode_)];
  travel_type_ = costing_->travel_type();
  access_mode_ = costing_->access_mode();

  // Initialize - create adjacency list, edgestatus support, A*, etc.
  Init(origin.edges.front().projected, destination.edges.front().projected);

  // Set origin and destination locations - seeds the adj. lists
  // Note: because we can correlate to more than one place for a given
  // PathLocation using edges.front here means we are only setting the
  // heuristics to one of them alternate paths using the other correlated
  // points to may be harder to find
  SetOrigin(graphreader, origin);
  SetDestination(graphreader, destination);

  // Find shortest path. Switch between a forward direction and a reverse
  // direction search based on the current costs. Alternating like this
  // prevents one tree from expanding much more quickly (if in a sparser
  // portion of the graph) rather than strictly alternating.
  // TODO - CostMatrix alternates, maybe should try alternating here?
  int n = 1;
  uint32_t forward_pred_idx, reverse_pred_idx;
  BDEdgeLabel pred, pred2;
  const GraphTile* tile;
  const GraphTile* tile2;
  bool expand_forward  = true;
  bool expand_reverse  = true;
  while (true) {
    // Allow this process to be aborted
    if (interrupt && (n % kInterruptIterationsInterval) == 0) {
      (*interrupt)();
    }
    n++;

    // Get the next predecessor (based on which direction was
    // expanded in prior step)
    if (expand_forward) {
      forward_pred_idx = forward_labels_.pop();
      if (forward_pred_idx != kInvalidLabel) {
        // Check if the edge on the forward search connects to a
        // settled edge on the reverse search tree.
        pred = forward_labels_.edgelabel(forward_pred_idx);
        if (reverse_labels_.edgestatus(pred.opp_edgeid()).set() == EdgeSet::kPermanent) {
          SetForwardConnection(pred);
        }
      } else {
        // Search is exhausted. If a connection has been found, return it
        if (best_connection_.cost < std::numeric_limits<float>::max()) {
          return FormPath(graphreader);
        } else {
          // No route found.
          LOG_ERROR("Bi-directional route failure - forward search exhausted: n = " +
                  std::to_string(forward_labels_.size()) + "," +
                  std::to_string(reverse_labels_.size()));
          return { };
        }
      }
    }
    if (expand_reverse) {
      reverse_pred_idx = reverse_labels_.pop();
      if (reverse_pred_idx != kInvalidLabel) {
        // Check if the edge on the reverse search connects to a
        // settled edge on the forward search tree.
        pred2 = reverse_labels_.edgelabel(reverse_pred_idx);
        if (forward_labels_.edgestatus(pred2.opp_edgeid()).set() == EdgeSet::kPermanent) {
          SetReverseConnection(pred2);
        }
      } else {
        // Search is exhausted. If a connection has been found, return it
        if (best_connection_.cost < std::numeric_limits<float>::max()) {
          return FormPath(graphreader);
        } else {
          // No route found.
          LOG_ERROR("Bi-directional route failure - reverse search exhausted: n = " +
                  std::to_string(forward_labels_.size()) + "," +
                  std::to_string(reverse_labels_.size()));
          return { };
        }
      }
    }

    // Terminate some number of iterations after an initial connection
    // has been found. This is not ideal, probably needs to be based on
    // the max edge cost but that has performance limitations,
    // so for now we use this bit of a hack...stay tuned.
    if (best_connection_.cost < std::numeric_limits<float>::max()) {
      if (forward_labels_.size() + reverse_labels_.size() > threshold_) {
        return FormPath(graphreader);
      }
    }

    // Expand from the search direction with lower sort cost.
    if ((pred.sortcost() + cost_diff_) < pred2.sortcost()) {
      // Expand forward - set to get next edge from forward adj. list
      // on the next pass
      expand_forward = true;
      expand_reverse = false;

      // Settle this edge.
      forward_labels_.settle(pred.edgeid());

      // Prune path if predecessor is not a through edge or if the maximum
      // number of upward transitions has been exceeded on this hierarchy level.
      if ((pred.not_thru() && pred.not_thru_pruning()) ||
          hierarchy_limits_forward_[pred.endnode().level()].StopExpanding()) {
        continue;
      }

      // Expand from the end node in forward direction.
      ExpandForward(graphreader, pred.endnode(), pred, forward_pred_idx, false);
    } else {
      // Expand reverse - set to get next edge from reverse adj. list
      // on the next pass
      expand_forward = false;
      expand_reverse = true;

      // Settle this edge
      reverse_labels_.settle(pred2.edgeid());

      // Prune path if predecessor is not a through edge
      if ((pred2.not_thru() && pred2.not_thru_pruning()) ||
          hierarchy_limits_reverse_[pred2.endnode().level()].StopExpanding()) {
        continue;
      }

      // Get the opposing predecessor directed edge. Need to make sure we get
      // the correct one if a transition occurred
      const DirectedEdge* opp_pred_edge =
        graphreader.GetGraphTile(pred2.opp_edgeid())->directededge(pred2.opp_edgeid());

      // Expand from the end node in reverse direction.
      ExpandReverse(graphreader, pred2.endnode(), pred2, reverse_pred_idx,
                    opp_pred_edge, false);
    }
  }
  return {};    // If we are here the route failed
}

// The edge on the forward search connects to a reached edge on the reverse
// search tree. Check if this is the best connection so far and set the
// search threshold.
void BidirectionalAStar::SetForwardConnection(const BDEdgeLabel& pred) {
  // Disallow connections that are part of a complex restriction.
  // TODO - validate that we do not need to "walk" the paths forward
  // and backward to see if they match a restriction.
  if (pred.on_complex_rest()) {
    return;
  }

  // Set a threshold to extend search
  if (threshold_ == 0) {
    threshold_ = GetThreshold(mode_, forward_labels_.size() + reverse_labels_.size());
  }

  GraphId oppedge = pred.opp_edgeid();
  EdgeStatusInfo oppedgestatus = reverse_labels_.edgestatus(oppedge);
  uint32_t predidx = reverse_labels_.edgelabel(oppedgestatus.index()).predecessor();
  float oppcost = (predidx == kInvalidLabel) ?
        0 : reverse_labels_.edgelabel(predidx).cost().cost;
  float c = pred.cost().cost + oppcost +
      reverse_labels_.edgelabel(oppedgestatus.index()).transition_cost();
  if (c < best_connection_.cost) {
    best_connection_ = { pred.edgeid(), oppedge, c };
  }
}

// The edge on the reverse search connects to a reached edge on the forward
// search tree. Check if this is the best connection so far and set the
// search threshold.
void BidirectionalAStar::SetReverseConnection(const BDEdgeLabel& pred) {
  // Disallow connections that are part of a cmplex restriction.
  // TODO - validate that we do not need to "walk" the paths forward
  // and backward to see if they match a restriction.
  if (pred.on_complex_rest()) {
    return;
  }

  // Set a threshold to extend search
  if (threshold_ == 0) {
    threshold_ = GetThreshold(mode_, forward_labels_.size() + reverse_labels_.size());
  }

  // Get the opposing edge - if this edge has been reached then a shortest
  // path has been found to the end node of this directed edge.
  GraphId oppedge = pred.opp_edgeid();
  EdgeStatusInfo oppedgestatus = forward_labels_.edgestatus(oppedge);
  uint32_t predidx = forward_labels_.edgelabel(oppedgestatus.index()).predecessor();
  float oppcost = (predidx == kInvalidLabel) ?
        0 : forward_labels_.edgelabel(predidx).cost().cost;
  float c = pred.cost().cost + oppcost +
      forward_labels_.edgelabel(oppedgestatus.index()).transition_cost();
  if (c < best_connection_.cost) {
    best_connection_ = { oppedge, pred.edgeid(), c };
  }
}

// Add edges at the origin to the forward adjacency list.
void BidirectionalAStar::SetOrigin(GraphReader& graphreader,
                 PathLocation& origin) {
  // Only skip inbound edges if we have other options
  bool has_other_edges = false;
  std::for_each(origin.edges.cbegin(), origin.edges.cend(), [&has_other_edges](const PathLocation::PathEdge& e){
    has_other_edges = has_other_edges || !e.end_node();
  });

  // Iterate through edges and add to adjacency list
  const NodeInfo* nodeinfo = nullptr;
  for (const auto& edge : origin.edges) {
    // If origin is at a node - skip any inbound edge (dist = 1)
    if (has_other_edges && edge.end_node()) {
      continue;
    }

    // Get the directed edge
    GraphId edgeid = edge.id;
    const GraphTile* tile = graphreader.GetGraphTile(edgeid);
    const DirectedEdge* directededge = tile->directededge(edgeid);

    // Get the tile at the end node. Skip if tile not found as we won't be
    // able to expand from this origin edge.
    const GraphTile* endtile = graphreader.GetGraphTile(directededge->endnode());
    if (endtile == nullptr) {
      continue;
    }

    // Get cost and sort cost (based on distance from endnode of this edge
    // to the destination
    nodeinfo = endtile->node(directededge->endnode());
    Cost cost = costing_->EdgeCost(directededge) * (1.0f - edge.dist);

    // We need to penalize this location based on its score (distance in meters from input)
    // We assume the slowest speed you could travel to cover that distance to start/end the route
    // TODO: assumes 1m/s which is a maximum penalty this could vary per costing model
    cost.cost += edge.score;
    float dist = astarheuristic_forward_.GetDistance(nodeinfo->latlng());
    float sortcost = cost.cost + astarheuristic_forward_.Get(dist);

    // Add edge to the LabelSet
    forward_labels_.add(edgeid, {kInvalidLabel, edgeid, directededge, cost,
            sortcost, dist, mode_});

    // Set the initial not_thru flag to false. There is an issue with not_thru
    // flags on small loops. Set this to false here to override this for now.
    forward_labels_.edgelabels().back().set_not_thru(false);
  }
}

// Add destination edges to the reverse path adjacency list.
void BidirectionalAStar::SetDestination(GraphReader& graphreader,
                     const PathLocation& dest) {
  // Only skip outbound edges if we have other options
  bool has_other_edges = false;
  std::for_each(dest.edges.cbegin(), dest.edges.cend(), [&has_other_edges](const PathLocation::PathEdge& e){
    has_other_edges = has_other_edges || !e.begin_node();
  });

  // Iterate through edges and add to adjacency list
  Cost c;
  for (const auto& edge : dest.edges) {
    // If the destination is at a node, skip any outbound edges (so any
    // opposing inbound edges are not considered)
    if (has_other_edges && edge.begin_node()) {
      continue;
    }

    // Get the directed edge
    GraphId edgeid = edge.id;
    const GraphTile* tile = graphreader.GetGraphTile(edgeid);
    const DirectedEdge* directededge = tile->directededge(edgeid);

    // Get the opposing directed edge, continue if we cannot get it
    GraphId opp_edge_id = graphreader.GetOpposingEdgeId(edgeid);
    if (!opp_edge_id.Is_Valid()) {
      continue;
    }
    const DirectedEdge* opp_dir_edge = graphreader.GetOpposingEdge(edgeid);

    // Get cost and sort cost (based on distance from endnode of this edge
    // to the origin. Make sure we use the reverse A* heuristic. Use the
    // directed edge for costing, as this is the forward direction along the
    // destination edge. Note that the end node of the opposing edge is in the
    // same tile as the directed edge.
    Cost cost = costing_->EdgeCost(directededge) * edge.dist;

    // We need to penalize this location based on its score (distance in meters from input)
    // We assume the slowest speed you could travel to cover that distance to start/end the route
    // TODO: assumes 1m/s which is a maximum penalty this could vary per costing model
    cost.cost += edge.score;
    float dist = astarheuristic_reverse_.GetDistance(tile->node(
                    opp_dir_edge->endnode())->latlng());
    float sortcost = cost.cost + astarheuristic_reverse_.Get(dist);

    // Add edge to the LabelSet
    reverse_labels_.add(opp_edge_id, {kInvalidLabel, opp_edge_id, edgeid,
        opp_dir_edge, cost, sortcost, dist, mode_, c, false});

    // Set the initial not_thru flag to false. There is an issue with not_thru
    // flags on small loops. Set this to false here to override this for now.
    reverse_labels_.edgelabels().back().set_not_thru(false);
  }
}

// Form the path from the adjacency list.
std::vector<PathInfo> BidirectionalAStar::FormPath(GraphReader& graphreader) {
  // Get the indexes where the connection occurs.
  uint32_t idx1 = forward_labels_.edgestatus(best_connection_.edgeid).index();
  uint32_t idx2 = reverse_labels_.edgestatus(best_connection_.opp_edgeid).index();

  // Metrics (TODO - more accurate cost)
  uint32_t pathcost = forward_labels_.edgelabel(idx1).cost().cost +
                      reverse_labels_.edgelabel(idx2).cost().cost;
  LOG_DEBUG("path_cost::" + std::to_string(pathcost));
  LOG_DEBUG("FormPath path_iterations::" + std::to_string(forward_labels_.size()) +
           "," + std::to_string(reverse_labels_.size()));

  // Work backwards on the forward path
  std::vector<PathInfo> path;
  for (auto edgelabel_index = idx1; edgelabel_index != kInvalidLabel;
      edgelabel_index = forward_labels_.edgelabel(edgelabel_index).predecessor()) {
    const BDEdgeLabel& edgelabel = forward_labels_.edgelabel(edgelabel_index);
    path.emplace_back(edgelabel.mode(), edgelabel.cost().secs,
                      edgelabel.edgeid(), 0);
  }

  // Reverse the list
  std:reverse(path.begin(), path.end());

  // Special case code if the last edge of the forward path is
  // the destination edge - update the elapsed time
  if (reverse_labels_.edgelabel(idx2).predecessor() == kInvalidLabel) {
    if (path.size() > 1) {
      path.back().elapsed_time = path[path.size()-2].elapsed_time +
          reverse_labels_.edgelabel(idx2).cost().secs;
    } else {
      path.back().elapsed_time = reverse_labels_.edgelabel(idx2).cost().secs;
    }
    return path;
  }

  // Get the elapsed time at the end of the forward path. NOTE: PathInfo
  // stores elapsed time as uint32_t but EdgeLabels uses float. Need to
  // accumulate in float and cast to int so we do not accumulate roundoff
  // error.
  float secs = path.back().elapsed_time;

  // Get the transition cost at the last edge of the reverse path
  float tc = reverse_labels_.edgelabel(idx2).transition_secs();

  // Append the reverse path from the destination - use opposing edges
  // The first edge on the reverse path is the same as the last on the forward
  // path, so get the predecessor.
  uint32_t edgelabel_index = reverse_labels_.edgelabel(idx2).predecessor();
  while (edgelabel_index != kInvalidLabel) {
    const BDEdgeLabel& edgelabel = reverse_labels_.edgelabel(edgelabel_index);
    GraphId oppedge = graphreader.GetOpposingEdgeId(edgelabel.edgeid());

    // Get elapsed time on the edge, then add the transition cost at
    // prior edge.
    uint32_t predidx = edgelabel.predecessor();
    if (predidx == kInvalidLabel) {
      secs += edgelabel.cost().secs;
    } else {
      secs += edgelabel.cost().secs - reverse_labels_.edgelabel(predidx).cost().secs;
    }
    secs += tc;
    path.emplace_back(edgelabel.mode(), secs, oppedge, 0);

    // Update edgelabel_index and transition cost to apply at next iteration
    edgelabel_index = predidx;
    tc = edgelabel.transition_secs();
  }
  return path;
}

}
}
