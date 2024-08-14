#ifndef FLSP_TOPO_UNSTRUCTURED_CONFIG_HH
#define FLSP_TOPO_UNSTRUCTURED_CONFIG_HH

#include <cstddef>

#include <flecsi/topo/unstructured/interface.hh>

namespace flsp::topo::unstructured {

template<std::size_t D>
struct config;

// 1D
template<>
struct config<1> : flecsi::topo::help {
  static constexpr std::size_t dimension() {
    return 1;
  }
  enum index_space {
    vertices,
    edges = vertices,
    faces = edges,
    cells,
    sides = cells,
    corners
  };
  using index_spaces = has<vertices, cells, corners>;
  using connectivities = list<from<cells, to<vertices, corners>>,
    from<vertices, to<cells, corners>>>;
  enum entity_list { owned, exclusive, shared, ghost };
  using entity_lists =
    list<entity<vertices, has<owned, exclusive, shared, ghost>>,
      entity<cells, has<owned, exclusive, shared, ghost>>,
      entity<corners, has<owned, exclusive, shared, ghost>>>;
};

// 2D
template<>
struct config<2> : flecsi::topo::help {
  static constexpr std::size_t dimension() {
    return 2;
  }
  enum index_space { vertices, edges, faces = edges, cells, sides, corners };
  using index_spaces = has<vertices, edges, cells, sides, corners>;
  using connectivities = list<from<cells, to<faces, vertices, sides, corners>>,
    from<vertices, to<cells, edges, sides, corners>>,
    from<faces, to<cells, vertices, sides>>,
    from<sides, to<cells, faces, vertices, corners>>,
    from<corners, to<cells, vertices, sides, edges>>>;
  enum entity_list { owned, exclusive, shared, ghost };
  using entity_lists =
    list<entity<vertices, has<owned, exclusive, shared, ghost>>,
      entity<cells, has<owned, exclusive, shared, ghost>>,
      entity<faces, has<owned, exclusive, shared, ghost>>,
      entity<sides, has<owned, exclusive, shared, ghost>>,
      entity<corners, has<owned, exclusive, shared, ghost>>>;
};

// 3D
template<>
struct config<3> : flecsi::topo::help {
  static constexpr std::size_t dimension() {
    return 3;
  }
  enum index_space { vertices, edges, faces, cells, sides, corners };
  using index_spaces = has<vertices, edges, faces, cells, sides, corners>;
  using connectivities =
    list<from<cells, to<edges, faces, vertices, sides, corners>>,
      from<vertices, to<cells, edges, sides, corners>>,
      from<edges, to<cells, vertices>>,
      from<faces, to<cells, vertices, sides>>,
      from<sides, to<cells, faces, vertices, corners>>,
      from<corners, to<cells, vertices, sides, edges>>>;
  enum entity_list { owned, exclusive, shared, ghost };
  using entity_lists =
    list<entity<vertices, has<owned, exclusive, shared, ghost>>,
      entity<cells, has<owned, exclusive, shared, ghost>>,
      entity<edges, has<owned, exclusive, shared, ghost>>,
      entity<faces, has<owned, exclusive, shared, ghost>>,
      entity<sides, has<owned, exclusive, shared, ghost>>,
      entity<corners, has<owned, exclusive, shared, ghost>>>;
};

template<std::size_t D>
using entity_kind = typename config<D>::index_space;

template<std::size_t D>
constexpr auto interface_kind();
template<>
constexpr auto
interface_kind<2>() {
  return entity_kind<2>::edges;
}
template<>
constexpr auto
interface_kind<3>() {
  return entity_kind<3>::faces;
}

template<std::size_t D, entity_kind<D> K>
constexpr auto
kind_id() -> std::enable_if_t<D == 2, std::uint32_t> {
  switch(K) {
    case entity_kind<D>::cells:
      return 0;
      break;
    case entity_kind<D>::vertices:
      return 1;
      break;
    case entity_kind<D>::edges:
      return 2;
      break;
    case entity_kind<D>::sides:
      return 3;
      break;
    case entity_kind<D>::corners:
      return 4;
      break;
  } // switch
} // kind_id

template<std::size_t D, entity_kind<D> K>
constexpr auto
kind_id() -> std::enable_if_t<D == 3, std::uint32_t> {
  switch(K) {
    case entity_kind<D>::cells:
      return 0;
      break;
    case entity_kind<D>::vertices:
      return 1;
      break;
    case entity_kind<D>::edges:
      return 2;
      break;
    case entity_kind<D>::faces:
      return 3;
      break;
    case entity_kind<D>::sides:
      return 4;
      break;
    case entity_kind<D>::corners:
      return 5;
      break;
  } // switch
} // kind_id

template<std::size_t D>
constexpr std::uint32_t num_index_spaces();
template<>
constexpr std::uint32_t
num_index_spaces<2>() {
  return 5;
};
template<>
constexpr std::uint32_t
num_index_spaces<3>() {
  return 6;
};

} // namespace flsp::topo::unstructured

#endif // FLSP_TOPO_UNSTRUCTURED_CONFIG_HH
