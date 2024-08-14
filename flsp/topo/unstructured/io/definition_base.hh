#ifndef FLSP_TOPO_UNSTRUCTURED_IO_DEFINITION_BASE_HH
#define FLSP_TOPO_UNSTRUCTURED_IO_DEFINITION_BASE_HH

#include "flsp/topo/unstructured/config.hh"
#include "flsp/topo/unstructured/util/common.hh"

#include <flecsi/topo/unstructured/types.hh>
#include <flecsi/util/array_ref.hh>
#include <flecsi/util/serialize.hh>

#include <optional>
#include <tuple>
#include <vector>

namespace flsp::topo::unstructured::io {

template<typename E, E... Ks>
struct required_keys {
  static constexpr bool value{true};
};

struct face_info {
  util::crs c2f;
  util::crs f2v;
  std::vector<util::gid> p2m;
  friend std::ostream & operator<<(std::ostream & stream, face_info const & d) {
    stream << "c2f:\n"
           << d.c2f << "\nf2v:\n"
           << d.f2v << "\np2m:\n"
           << flecsi::flog::container{d.p2m} << std::endl;
    return stream;
  } // operator<<
};
// Use this to see ids without ones' complement.
struct debug_face_info : face_info {
  debug_face_info(face_info const & fi) : face_info(fi) {}
  friend std::ostream & operator<<(std::ostream & stream,
    debug_face_info const & d) {
    stream << "c2f:\n"
           // scrs removes ones' complement
           << util::scrs{d.c2f} << "\nf2v:\n"
           << d.f2v << "\np2m:\n"
           << flecsi::flog::container{d.p2m} << std::endl;
    return stream;
  } // operator<<
};

using mat_ids = std::vector<std::vector<std::pair<util::gid, float>>>;
using bnd_ids = std::vector<std::vector<util::gid>>;

using iota_view = flecsi::util::iota_view<std::size_t>;

/*!
  Abstract base class for mesh definitions.
 */

template<std::size_t D>
struct definition_base {

  virtual ~definition_base(){};

  /*!
    Return the global number of entities of the given kind.
   */

  virtual util::gid num_entities(entity_kind<D> k) const = 0;

  /*!
    Return relational information for the given cell range.

    @param r The range of cells expressed as a FleCSI iota_view.

    @return A std::tuple containing cell-to-vertex connectivity (crs),
      face connectivities (face_info), and materials IDs (mat_ids).
   */

  virtual std::
    tuple<util::crs, std::optional<face_info>, std::optional<mat_ids>>
    cell_data(iota_view const & r) const = 0;

  /*!
    Return the requested range of coordinates.

    @param r The range of vertices expressed as a FleCSI iota_view.

    @return A std::tuple containing coordinates (std::vector<point>)
      and boundary IDs (bnd_ids).
   */

  virtual std::tuple<std::vector<util::point<D>>, std::optional<bnd_ids>>
  vertex_data(iota_view const & r) const = 0;

}; // definition_base

/*!
  Empty mesh definition class to use for non-root process construction.
 */

template<std::size_t D>
struct undefined_definition : definition_base<D> {
  undefined_definition() {}
  util::gid num_entities(entity_kind<D>) const override {
    flog_fatal("undefined mesh definition");
    return {};
  }
  std::tuple<util::crs, std::optional<face_info>, std::optional<mat_ids>>
  cell_data(iota_view const &) const override {
    flog_fatal("undefined mesh definition");
    return {};
  }
  std::tuple<std::vector<util::point<D>>, std::optional<bnd_ids>> vertex_data(
    iota_view const &) const override {
    flog_fatal("undefined mesh definition");
    return {};
  }
}; // struct undefined_definition

} // namespace flsp::topo::unstructured::io

namespace flecsi::util::serial {
template<class T>
struct traits<std::optional<T>, std::enable_if_t<!bit_copyable_v<T>>> {
  using type = std::optional<T>;
  template<class P>
  static void put(P & p, const type & v) {
    if(v.has_value()) {
      serial::put(p, bool(true));
      serial::put<T>(p, v.value());
    }
    else {
      serial::put(p, bool(false));
    }
  }
  static type get(const std::byte *& p) {
    bool has_value = serial::get<bool>(p);
    if(has_value) {
      return std::make_optional(std::move(serial::get<T>(p)));
    }
    return std::nullopt;
  }
};

template<>
struct traits<flsp::topo::unstructured::io::face_info> {
  using type = flsp::topo::unstructured::io::face_info;
  template<class P>
  static void put(P & p, const type & fi) {
    serial::put(p, fi.c2f, fi.f2v, fi.p2m);
  }
  static type get(const std::byte *& p) {
    const cast r{p};
    return type{r, r, r};
  }
};
} // namespace flecsi::util::serial

#endif // FLSP_TOPO_UNSTRUCTURED_IO_DEFINITION_BASE_HH
