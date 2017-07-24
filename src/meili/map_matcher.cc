#include <cmath>
#include "midgard/distanceapproximator.h"
#include "loki/search.h"
#include "meili/routing.h"
#include "meili/geometry_helpers.h"
#include "meili/match_route.h"
#include "meili/map_matcher.h"


namespace {

using namespace valhalla;
using namespace valhalla::meili;


inline float
GreatCircleDistanceSquared(const Measurement& left,
                           const Measurement& right)
{ return left.lnglat().DistanceSquared(right.lnglat()); }


inline float
GreatCircleDistance(const Measurement& left,
                    const Measurement& right)
{ return left.lnglat().Distance(right.lnglat()); }

inline float
ClockDistance(const Measurement& left,
              const Measurement& right)
{ return right.epoch_time() - left.epoch_time(); }


struct Interpolation {
  midgard::PointLL projected;
  baldr::GraphId edgeid;
  float sq_distance;
  float route_distance;
  float route_time;
  float edge_distance;

  float sortcost(const MapMatching& mm, float gc_dist, float clk_dist) const
  { return mm.CalculateTransitionCost(0.f, route_distance, gc_dist, route_time, clk_dist) +
      mm.CalculateEmissionCost(sq_distance); }
};


// Find the interpolation along the route where the transition cost +
// emission cost is minimal
template <typename segment_iterator_t>
Interpolation
InterpolateMeasurement(const MapMatching& mapmatching,
                       segment_iterator_t begin,
                       segment_iterator_t end,
                       const Measurement& measurement,
                       float match_measurement_distance,
                       float match_measurement_time)
{
  const baldr::GraphTile* tile(nullptr);
  midgard::DistanceApproximator approximator(measurement.lnglat());

  // Route distance from each segment begin to the beginning segment
  float segment_begin_route_distance = 0.f;
  float minimal_cost = std::numeric_limits<float>::infinity();

  // Invalid edgeid indicates that no interpolation found
  Interpolation best_interp;

  for (auto segment = begin; segment != end; segment++) {
    const auto directededge = mapmatching.graphreader().directededge(segment->edgeid, tile);
    if (!directededge) {
      continue;
    }

    const auto edgeinfo = tile->edgeinfo(directededge->edgeinfo_offset());

    const auto& shape = edgeinfo.shape();
    if (shape.empty()) {
      continue;
    }

    midgard::PointLL projected_point;
    float sq_distance, offset;
    std::tie(projected_point, sq_distance, std::ignore, offset) = helpers::Project(measurement.lnglat(), shape, approximator);

    // Find out the correct offset
    if (!directededge->forward()) {
      offset = 1.f - offset;
    }

    // Distance from the projected point to the segment begin, or
    // segment begin if we walk reversely
    const auto distance_to_segment_ends = std::abs(directededge->length() * (offset - segment->source));

    // The absolute route distance from projected point to the
    // beginning segment
    const auto route_distance = segment_begin_route_distance + distance_to_segment_ends;

    // Get the amount of time spent on this segment
    auto edge_percent = segment->target - segment->source;
    auto route_time = mapmatching.costing()->EdgeCost(directededge).secs * edge_percent;

    Interpolation interp{projected_point, segment->edgeid, sq_distance, route_distance, route_time, offset};
    const auto cost = interp.sortcost(mapmatching, match_measurement_distance, match_measurement_time);
    if (cost < minimal_cost) {
      minimal_cost = cost;
      best_interp = std::move(interp);
    }

    // Assume segments are connected
    segment_begin_route_distance += directededge->length() * edge_percent;
  }

  return best_interp;
}


// Interpolate measurements along the route from previous state to
// current state, or along the route from current state to next state
std::vector<MatchResult>
InterpolateMeasurements(const MapMatching& mapmatching,
                        const MapMatching::state_iterator& previous_state,
                        const MapMatching::state_iterator& state,
                        const MapMatching::state_iterator& next_state,
                        const std::vector<Measurement>& interpolated_measurements)
{
  //nothing to do here
  if (interpolated_measurements.empty())
    return {};

  //can't interpolate these because they don't happen between two valid states or
  //we weren't able to get a downstream route
  std::vector<MatchResult> results;
  std::vector<EdgeSegment> route;
  if (!state.IsValid() || !next_state.IsValid() ||
    MergeRoute(route, *state, *next_state).empty()) {
    for (const auto& measurement: interpolated_measurements)
      results.emplace_back(MatchResult{measurement.lnglat(), 0.f, baldr::GraphId{}, -1.f, measurement.epoch_time(), kInvalidStateId});
    return results;
  }

  //for each point that needs interpolated
  for (const auto& measurement: interpolated_measurements) {
    const auto& match_measurement = mapmatching.measurement(state->time());
    const auto match_measurement_distance = GreatCircleDistance(measurement, match_measurement);
    const auto match_measurement_time = ClockDistance(measurement, match_measurement);
    //interpolate this point along the route
    const auto& interp = InterpolateMeasurement(mapmatching, route.begin(), route.end(),
        measurement, match_measurement_distance, match_measurement_time);

    //if it was able to do the interpolation
    if (interp.edgeid.Is_Valid()) {
      //dont allow subsequent points to get interpolated before this point
      //we do this by editing the route to start where this point was interpolated
      auto itr = std::find_if(route.begin(), route.end(), [&interp](const EdgeSegment& e){ return e.edgeid == interp.edgeid; });
      itr = std::find_if(itr, route.end(), [&interp](const EdgeSegment& e){ return e.edgeid == interp.edgeid && e.target > interp.edge_distance; });
      if(itr != route.end()) {
        itr->source = interp.edge_distance;
        route.erase(route.begin(), itr);
      }
      //keep the interpolated match result
      results.emplace_back(MatchResult{interp.projected, std::sqrt(interp.sq_distance), interp.edgeid, interp.edge_distance, measurement.epoch_time(), kInvalidStateId});
    }//couldnt interpolate this point
    else
      results.emplace_back(MatchResult{measurement.lnglat(), 0.f, baldr::GraphId{}, -1.f, measurement.epoch_time(), kInvalidStateId});
  }

  return results;
}


// Interplolate measurements grouped by previous match measurement
// time
std::unordered_map<Time, std::vector<MatchResult>>
InterpolateTimedMeasurements(const MapMatching& mapmatching,
                             const MapMatching::state_iterator& state_rbegin,
                             const MapMatching::state_iterator& state_rend,
                             const std::unordered_map<Time, std::vector<Measurement>>& interpolated_measurements,
                             Time time)
{
  std::unordered_map<Time, std::vector<MatchResult>> resultmap;

  for (auto previous_state = std::next(state_rbegin),
                     state = state_rbegin,
                next_state = MapMatching::state_iterator(nullptr);
       state != state_rend;
       next_state = state, state++, previous_state = std::next(state)) {

    const auto it = interpolated_measurements.find(time);
    if (it != interpolated_measurements.end()) {
      if (state.IsValid()) {
        const auto& results = InterpolateMeasurements(mapmatching, previous_state, state, next_state, it->second);
        resultmap.emplace(time, results);
      } else {
        auto& results = resultmap[time];
        for (const auto& measurement: it->second) {
          results.push_back({measurement.lnglat(), 0.f, {}, -1.f, measurement.epoch_time(), kInvalidStateId});
        }
      }
    }

    time--;
  }

  return resultmap;
}


// Find the match result of a state, given its previous state and next
// state
MatchResult
FindMatchResult(const MapMatching::state_iterator& previous_state,
                const State& state, const Measurement& measurement,
                const MapMatching::state_iterator& next_state)
{
  baldr::GraphId edgeid;

  // Construct the route from previous state to current state, and
  // find out which edge it matches
  if (previous_state.IsValid()) {
    // It must stay on the last edge of the route
    const auto rbegin = previous_state->RouteBegin(state),
                 rend = previous_state->RouteEnd();
    if (rbegin != rend) {
      edgeid = rbegin->edgeid;
    }
  }

  // Do the same from current state to next state
  if (!edgeid.Is_Valid() && next_state.IsValid()) {
    // It must stay on the first edge of the route
    for (auto label = state.RouteBegin(*next_state);
         label != state.RouteEnd(); label++) {
      if (label->edgeid.Is_Valid()) {
        edgeid = label->edgeid;
      }
    }
  }

  // If we get a valid edge and find it in the state that's good
  for(const auto& edge : state.candidate().edges)
    if(edge.id == edgeid)
      return {edge.projected, std::sqrt(edge.score), edgeid, edge.dist, measurement.epoch_time(), state.id()};

  // If we failed to get a valid edge or can't find it
  // At least we know which point it matches
  const auto& edge = state.candidate().edges.front();
  return {edge.projected, std::sqrt(edge.score), edgeid, edge.dist, measurement.epoch_time(), state.id()};
}


// Find the corresponding match results of a list of states
std::vector<MatchResult>
FindMatchResults(const MapMatching& mapmatching,
                 const MapMatching::state_iterator& state_rbegin,
                 const MapMatching::state_iterator& state_rend,
                 Time time)
{
  if (state_rbegin == state_rend) {
    return {};
  }

  std::vector<MatchResult> results;

  for (auto previous_state = std::next(state_rbegin),
                     state = state_rbegin,
                next_state = MapMatching::state_iterator(nullptr);
       state != state_rend;
       next_state = state, state++, previous_state = std::next(state)) {
    const auto& measurement = mapmatching.measurement(time);
    if (state.IsValid()) {
      results.push_back(FindMatchResult(previous_state, *state, measurement, next_state));
    } else {
      results.emplace_back(MatchResult{measurement.lnglat(), 0.f, baldr::GraphId{}, -1.f, measurement.epoch_time(), kInvalidStateId});
    }
    time--;
  }

  std::reverse(results.begin(), results.end());
  return results;
}


// Insert the interpolated results into the result list, and
// return a new result list
std::vector<MatchResult>
MergeMatchResults(const std::vector<MatchResult>& results,
                  const std::unordered_map<Time, std::vector<MatchResult>>& interpolated_results)
{
  std::vector<MatchResult> merged_results;

  Time time = 0;
  for (const auto& result: results) {
    merged_results.push_back(result);

    const auto it = interpolated_results.find(time);
    if (it != interpolated_results.end()) {
      for (const auto& result: it->second) {
        merged_results.push_back(result);
      }
    }

    time++;
  }

  return merged_results;
}

}


namespace valhalla {
namespace meili {

MapMatcher::MapMatcher(const boost::property_tree::ptree& config,
                       baldr::GraphReader& graphreader,
                       CandidateQuery& candidatequery,
                       const sif::cost_ptr_t* mode_costing,
                       sif::TravelMode travelmode)
    : config_(config),
      graphreader_(graphreader),
      candidatequery_(candidatequery),
      mode_costing_(mode_costing),
      travelmode_(travelmode),
      mapmatching_(graphreader_, mode_costing_, travelmode_, config_),
      interrupt_(nullptr) {}


MapMatcher::~MapMatcher() {}


std::vector<MatchResult>
MapMatcher::OfflineMatch(
    const std::vector<Measurement>& measurements, uint32_t k) {
  mapmatching_.Clear();

  const auto begin = measurements.begin(),
               end = measurements.end();

  if (begin == end) {
    return {};
  }

  const auto interpolation_distance = config_.get<float>("interpolation_distance");
  const auto sq_interpolation_distance = interpolation_distance * interpolation_distance;
  std::unordered_map<Time, std::vector<Measurement>> interpolated_measurements;
  /*
  // Make the locations to pass to loki::search
  std::vector<baldr::Location> locations;
  for(const auto& measurement : measurements) {
    //first or last one are states for sure, also any that are far enough apart from previous
    if(&measurement == &measurements.front() || &measurement == &measurements.back() ||
       sq_interpolation_distance < locations.back().latlng_.DistanceSquared(measurement.lnglat())) {
      locations.emplace_back(baldr::Location{measurement.lnglat(), baldr::Location::StopType::BREAK, 0,
          static_cast<unsigned int>(std::max(measurement.search_radius(), measurement.gps_accuracy()))});
      //TODO: supplant measurement with a location object
      //for now use a string trick to store an index of the original measurement
      uint32_t index = locations.size() - 1;
      locations.back().city_.resize(4);
      memcpy(&locations.back().city_[0], static_cast<const char*>(static_cast<const void*>(&index)), 4);
    }//otherwise it gets interpolated
    else {
      interpolated_measurements[locations.size() - 1].push_back(measurement);
    }
  }

  // Get the candidates for each state
  auto candidates = loki::Search(locations, graphreader_,
      mapmatching_.costing()->GetEdgeFilter(), mapmatching_.costing()->GetNodeFilter());
  std::vector<baldr::PathLocation> correlated; correlated.reserve(candidates.size());
  for(const auto& candidate : candidates)
    correlated.push_back(candidate.second);
  std::sort(correlated.begin(), correlated.end(), [](const baldr::PathLocation& a, const baldr::PathLocation& b) {
    //TODO: stop abusing the location stuff to get an index out, see above
    uint32_t ai = *static_cast<const uint32_t*>(static_cast<const void*>(&a.city_[0]));
    uint32_t bi = *static_cast<const uint32_t*>(static_cast<const void*>(&b.city_[0]));
    return ai < bi;
  });

  for(const auto& c : correlated) {
    std::cout << "-------------" << std::endl;
    for(const auto& e : c.edges) {
      std::cout << e.id << "  " << e.dist << "  " << e.score/10.f << std::endl;
    }
    std::cout << "*************" << std::endl;
    const auto& candidates = candidatequery_.Query(c.latlng_,c.radius_ * c.radius_,mapmatching_.costing()->GetEdgeFilter());
    for(const auto& b : candidates) {
      for(const auto& e : b.edges) {
        std::cout << e.id << "  " << e.dist << "  " << e.score << std::endl;
      }
    }
    std::cout << "-------------" << std::endl;
  }


  // Add the states to the HMM
  for(auto& candidate : candidates) {
    //TODO: stop abusing the location stuff to get an index out, see above
    uint32_t index = *static_cast<const uint32_t*>(static_cast<const void*>(&candidate.first.city_[0]));
    const auto& measurement = measurements[index];
    //TODO: loki is hacked to add a multiplier to the score we undo that
    //here so that meili can rely on it being sq dist
    for(auto& edge : candidate.second.edges)
      edge.score /= 10.f;
    //actually add the state
    mapmatching_.AppendState(measurement, candidate.second);
  }
*/

  // Always match the first measurement
  auto time = AppendMeasurement(*begin);
  auto latest_match_measurement = begin;
  std::size_t match_count = 1;

  for (auto measurement = next(begin); measurement != end; measurement++) {
    const auto sq_distance = GreatCircleDistanceSquared(*latest_match_measurement, *measurement);
    // Always match the last measurement
    if (sq_interpolation_distance < sq_distance || std::next(measurement) == end) {
      time = AppendMeasurement(*measurement);
      latest_match_measurement = measurement;
      match_count++;
    } else {
      interpolated_measurements[time].push_back(*measurement);
    }
  }

  //const auto time = locations.size() - 1;
  //const auto match_count = locations.size();
  const auto state_rbegin = mapmatching_.SearchPath(time);
  const auto state_rend = std::next(state_rbegin, match_count);
  const auto& results = FindMatchResults(mapmatching_, state_rbegin, state_rend, time);

  // Done if no measurements to interpolate
  if (interpolated_measurements.empty()) {
    return results;
  }

  // Find results for all interpolated measurements
  const auto& interpolated_results = InterpolateTimedMeasurements(
      mapmatching_,
      state_rbegin,
      state_rend,
      interpolated_measurements,
      time);

  // Insert the interpolated results into the result list
  return MergeMatchResults(results, interpolated_results);
}


Time
MapMatcher::AppendMeasurement(const Measurement& measurement)
{
  // Test interrupt
  if (interrupt_) {
    (*interrupt_)();
  }
  const auto& candidates = candidatequery_.Query(
      measurement.lnglat(),
      std::max(measurement.sq_search_radius(), measurement.sq_gps_accuracy()),
      mapmatching_.costing()->GetEdgeFilter());
  return mapmatching_.AppendState(measurement, candidates.begin(), candidates.end());
  /*baldr::Location loc(measurement.lnglat(), baldr::Location::StopType::BREAK, 0,
      static_cast<unsigned int>(std::max(measurement.search_radius(), measurement.gps_accuracy())));
  auto candidate = loki::Search({loc}, graphreader_,
      mapmatching_.costing()->GetEdgeFilter(), mapmatching_.costing()->GetNodeFilter()).begin()->second;
  for(auto& e : candidate.edges)
    e.score /= 10.f;
  return mapmatching_.AppendState(measurement, candidate);*/
}

}
}
