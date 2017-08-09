#ifndef VALHALLA_THOR_LABELSET_H_
#define VALHALLA_THOR_LABELSET_H_

#include <unordered_map>
#include <valhalla/baldr/graphid.h>
#include <valhalla/baldr/double_bucket_queue.h>
#include <valhalla/sif/edgelabel.h>
#include <valhalla/thor/edgestatus.h>

namespace valhalla {
namespace thor {

/**
 * Edge Label set. Includes EdgeLabels, adjacency list, and edge status.
 * The type of EdgeLabel is set by the template instantiation. Adjacency
 * list (priority queue) uses a double bucket queue for fast sorting.
 */
template <class label_t>
class LabelSet {
 public:

  /**
   * Clear the label set information.
   */
  void clear() {
    edgelabels_.clear();
    adjacencylist_.reset();
    edgestatus_.reset();
  }

  /**
   * Initialize the LabelSet.
   * @param edgelabel_count  Edge label count to reserve.
   * @param bucket_size      Bucket size for double bucket queue.
   * @param bucket_count     Number of buckets in double bucket queue.
   * @param mincost          Minimum cost in the adjacency list.
   */
  void Init(const uint32_t edgelabel_count, const uint32_t bucket_size,
            const uint32_t bucket_count, const float mincost) {
    // Reserve size for edge labels - do this here rather than in constructor
    // to limit how much extra memory is used for persistent objects
    edgelabels_.reserve(edgelabel_count);

    // Set up lambdas to get sort costs
    const auto edgecost = [this](const uint32_t label) {
      return edgelabels_[label].sortcost();
    };

    // Construct adjacency list and edge status
    adjacencylist_.reset(new baldr::DoubleBucketQueue(mincost,
          bucket_count * bucket_size, bucket_size, edgecost));
    edgestatus_.reset(new EdgeStatus(edgelabel_count));
  }

  /**
   * Gets the next edge label (by index) from the adjacency list.
   * @return  Returns the index of the next EdgeLabel. Returns
   *          kInvalidLabel if the queue is empty.
   */
  uint32_t pop() {
    return adjacencylist_->pop();
  }

  /**
   * Returns the size of the EdgeLabel list.
   * @return Returns the size of the label list.
   */
  uint32_t size() const {
    return edgelabels_.size();
  }

  /**
   * Get the edge status.
   * @param  edgeid  GraphId of the edge.
   * @return Returns the edge status information (index and set).
   */
  EdgeStatusInfo edgestatus(const baldr::GraphId& edgeid) const {
    return edgestatus_->Get(edgeid);
  }

  /**
   * Return a const reference to the EdgeLabel at the specified index.
   * @param  index  Index of the edge label.
   * @return  Returns a const reference to the edge label.
   */
  const label_t& edgelabel(const uint32_t index) const {
    return edgelabels_[index];
  }

  /**
   * Returns a reference to the EdgeLabels vector.
   * @return  Returns a reference to the edge labels vector.
   */
  std::vector<label_t>& edgelabels() {
    return edgelabels_;
  }

  /**
   * Settle this edge. Mark as permanently labeled, meaning the lowest cost
   * has been found to the end of this edge.
   * @param edgeid  GraphId of the edge to settle.
   */
  void settle(const baldr::GraphId& edgeid) {
    edgestatus_->Update(edgeid, EdgeSet::kPermanent);
  }

  /**
   * Possibly lower the cost to the end of an edge already in the adjacency
   * list. Checks if cost is lower, if it is the adjacency list is updated
   * (cost in the priority queue is updated) and the label is updated with
   * new cost and predecessor information.
   * @param  index     Edge Label index.
   * @param  pred_idx  Predecessor edge index in EdgeLabels.
   * @param  newcost   Updated cost.
   * @param  tc        Updated turn cost.
   */
  void decrease(const uint32_t index, const uint32_t pred_idx,
                const sif::Cost& newcost, const sif::Cost& tc) {
    auto& lab = edgelabels_[index];
    if (newcost.cost <  lab.cost().cost) {
      float newsortcost = lab.sortcost() - (lab.cost().cost - newcost.cost);
      adjacencylist_->decrease(index, newsortcost);
      lab.Update(pred_idx, newcost, newsortcost, tc);
    }
  }

  /**
   * Add and edge label to the EdgeLabel vector and to adjacency list. Set the
   * edge status to temporary (and set its index).
   * @param  edgeid  GraphId of the edge to add.
   * @param  label   EdgeLabel to add.
   */
  void add(const baldr::GraphId& edgeid, label_t&& label) {
    uint32_t idx = edgelabels_.size();
    edgestatus_->Set(edgeid, EdgeSet::kTemporary, idx);
    edgelabels_.emplace_back(label);
    adjacencylist_->add(idx);
  }

 private:
  // Edge status. Mark edges that are in adjacency list or settled.
  std::shared_ptr<EdgeStatus> edgestatus_;

  // Vector of edge labels (requires access by index).
  std::vector<label_t> edgelabels_;

  // Adjacency list - approximate double bucket sort
  std::shared_ptr<baldr::DoubleBucketQueue> adjacencylist_;
};

}
}

#endif  // VALHALLA_THOR_EDGESTATUS_H_
