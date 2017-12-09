#ifndef VALHALLA_MEILI_STATE_H_
#define VALHALLA_MEILI_STATE_H_

#include <cstdint>
#include <functional>
#include <limits>

namespace valhalla {
namespace meili {

constexpr uint32_t kInvalidTime = std::numeric_limits<uint32_t>::max();

struct StateId
{
 public:
  using Time = uint32_t;

  using Id = uint32_t;

  StateId()
      : time_(kInvalidTime),
        id_(0) {
  }

  explicit StateId(const Time& time, Id id)
      : time_(time),
        id_(id) {
  }

  Time time() const
  { return time_; }

  Id id() const
  { return id_; }

  bool IsValid() const
  { return time_ != kInvalidTime; }

  bool operator ==(const StateId& rhs) const
  { return time() == rhs.time() && id() == rhs.id(); }

  bool operator !=(const StateId& rhs) const
  { return time() != rhs.time() || id() != rhs.id(); }

  // Create a single 64 bit value
  uint64_t value() const
  { return time_ | id_ << 32; }

 private:
  uint64_t time_ : 32;
  uint64_t id_   : 32;
};

}
}

// Extend the standard namespace to know how to hash StateIds
namespace std {
template <>
struct hash<valhalla::meili::StateId> {
  inline std::size_t operator()(const valhalla::meili::StateId& stateid) const {
    return static_cast<size_t>(stateid.value());
  }
};

template <>
struct hash<std::pair<valhalla::meili::StateId, valhalla::meili::StateId> > {
  inline std::size_t operator()(const std::pair<valhalla::meili::StateId, valhalla::meili::StateId>& couple) const {
    auto seed = static_cast<size_t>(couple.first.value());
    return static_cast<size_t>(couple.second.value()) + 0x9e3779b9 + (seed<<6) + (seed>>2);
  }
};
}

#endif // VALHALLA_MEILI_STATE_H_
