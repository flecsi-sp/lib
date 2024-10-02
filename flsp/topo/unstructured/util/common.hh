#ifndef FLSP_TOPO_UNSTRUCTURED_UTIL_COMMON_HH
#define FLSP_TOPO_UNSTRUCTURED_UTIL_COMMON_HH

#include <flecsi/util/color_map.hh>
#include <flecsi/util/common.hh>
#include <flecsi/util/crs.hh>
#include <flecsi/util/mpi.hh>

#include <cmath>

/// @cond core
namespace flsp::topo::unstructured {

// Import types from flecsi
using flecsi::Color;
namespace flog = flecsi::flog;

namespace util {

/// @addtogroup util
/// \{

// Import types from flecsi::util
using id = flecsi::util::id;
using gid = flecsi::util::gid;
using crs = flecsi::util::crs;
using equal_map = flecsi::util::equal_map;
namespace mpi = flecsi::util::mpi;
using flecsi::util::force_unique;
using flecsi::util::to_vector;
using flecsi::util::unique_each;

// Define standard typedefs
template<std::size_t D>
using point = std::array<double, D>;

/*!
  Test if the sign bit is set for ids stored using the ones' complement
  strategy for orientation information. When set, this implies that the entity
  referenced by @em oid has a negative orientation with respect to its
  connected entity type.

  @param oid The id to test.
  @return true if sign bit is set, false otherwise.
 */
template<typename T>
inline auto
sign_bit(T oid) {
  return std::signbit(static_cast<std::make_signed_t<T>>(oid));
}

/*!
  Return the raw id referenced by @em oid for ids stored using the ones'
  complement strategy for orientation information.

  @param oid The ones' complement id.
 */
template<typename T>
inline auto
get_id(T oid) {
  return sign_bit(oid) ? ~oid : oid;
}

/*!
  Return a std::pair<bool, T> identifying the orientation and raw id,
  respectively, for ids stored using the ones' complement strategy for
  orientation information.

  @code
  for(auto c: m.cells()) {
    for(auto f: m.faces(c)) {
      [sgn, id] = get_sign_id(f);
      q += sgn*flux[id];
    }
  }
  @endcode
 */
template<typename T>
inline auto
get_sign_id(T oid) {
  return std::make_pair(-2 * sign_bit(oid) + 1, get_id(oid));
}

/*!
  Wrapper type to strip ones' complement from ids of a crs instance.
 */
struct scrs : crs {
  scrs(crs const & g) : crs(g) {}
};

/*!
  Expand crs graph while stripping ones' complement from ids.
 */
inline std::string
expand(scrs const & graph) {
  std::stringstream stream;
  std::size_t r{0};
  for(const auto row : graph) {
    stream << r++ << ": <";
    bool first = true;
    for(const std::size_t i : row) {
      if(first)
        first = false;
      else
        stream << ",";
      stream << get_id(i);
    }
    stream << ">" << std::endl;
  }
  return stream.str();
}

/*!
  Output crs graph while stripping ones' complement from ids.
 */
inline std::ostream &
operator<<(std::ostream & stream, scrs const & graph) {
  stream << "crs offsets: ";
  for(auto o : graph.offsets.ends())
    stream << o << " ";
  stream << "\n\ncrs indices: ";
  for(auto o : graph.values) {
    stream << get_id(o) << " ";
  }
  stream << "\n\ncrs expansion:\n" << expand(graph) << std::endl;
  return stream << std::endl;
} // operator<<

/// \}
} // namespace util
} // namespace flsp::topo::unstructured

/// @endcond

#endif // FLSP_TOPO_UNSTRUCTURED_UTIL_COMMON_HH
