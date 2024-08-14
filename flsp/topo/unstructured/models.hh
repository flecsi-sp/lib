/*
  Copyright (c) 2022, Triad National Security, LLC.
  All rights reserved.
 */

#ifndef FLSP_TOPO_UNSTRUCTURED_MODELS_HH
#define FLSP_TOPO_UNSTRUCTURED_MODELS_HH

#include "flsp/topo/unstructured/config.hh"
#include "flsp/topo/unstructured/util/common.hh"

#include <cstdint>
#include <type_traits>

namespace flsp::topo::unstructured::io {

template<std::size_t D, entity_kind<D> K>
struct model {};

// Corner creation is defined with 2D kinds (same for 2D and 3D).

// clang-format off
struct tri3 {
  static constexpr int num_edges{3};
  static constexpr std::uint32_t edges[num_edges][2] = {
    {0, 1}, /* v0, v1 */
    {1, 2},
    {2, 0}
  };
  // Sides use edge definitions.
}; // struct tri3

struct quad4 {
  static constexpr int num_edges{4};
  static constexpr std::uint32_t edges[num_edges][2] = {
    {0, 1}, /* v0, v1 */
    {1, 2},
    {2, 3},
    {3, 0}
  };
  // Sides use edge definitions.
}; // struct quad4
// clang-format on

// 2D Edges.
template<std::size_t D, entity_kind<D> K>
auto
create_cell_entities(std::tuple<util::gid, util::id, util::id> const & cid,
  util::crs const & c2v,
  std::map<entity_kind<D>, util::crs> const &)
  -> std::enable_if_t<D == 2 && K == entity_kind<D>::edges,
    std::pair<util::crs, std::map<entity_kind<D>, util::crs>>> {
  auto [gid, id, cfaid] = cid;
  auto const & vertices = c2v[id];
  util::crs entities;
  std::map<entity_kind<D>, util::crs> a2a;

  switch(vertices.size()) {
    case 3:
      for(int e{0}; e < tri3::num_edges; ++e) {
        entities.add_row(
          {vertices[tri3::edges[e][0]], vertices[tri3::edges[e][1]]});
      } // for
      break;
    case 4:
      for(int e{0}; e < quad4::num_edges; ++e) {
        entities.add_row(
          {vertices[quad4::edges[e][0]], vertices[quad4::edges[e][1]]});
      } // for
      break;
    default:
      flog_fatal("invalid edge definition");
  } // switch

  return std::make_pair(std::move(entities), std::move(a2a));
} // create_cell_entities

// 2D Sides.
template<std::size_t D, entity_kind<D> K>
auto
create_cell_entities(std::tuple<util::gid, util::id, util::id> const & cid,
  util::crs const & c2v,
  std::map<entity_kind<D>, util::crs> const & aux)
  -> std::enable_if_t<D == 2 && K == entity_kind<D>::sides,
    std::pair<util::crs, std::map<entity_kind<D>, util::crs>>> {
  auto [gid, id, cfaid] = cid;
  auto const & vertices = c2v[id];
  util::crs entities;
  std::map<entity_kind<D>, util::crs> a2a;

  auto const & corners =
    aux.count(entity_kind<D>::corners)
      ? std::make_optional(aux.at(entity_kind<D>::corners)[cfaid])
      : std::nullopt;

  auto const & interfaces = aux.at(interface_kind<D>())[cfaid];

  switch(vertices.size()) {
    case 3:
      for(int e{0}; e < tri3::num_edges; ++e) {
        entities.add_row(
          {vertices[tri3::edges[e][0]], vertices[tri3::edges[e][1]], gid});

        if(corners.has_value()) {
          a2a[entity_kind<D>::corners].add_row(
            {(*corners)[tri3::edges[e][0]], (*corners)[tri3::edges[e][1]]});
        } // if
        a2a[interface_kind<D>()].add_row({interfaces[e]});
      } // for
      break;
    case 4:
      for(int e{0}; e < quad4::num_edges; ++e) {
        entities.add_row(
          {vertices[quad4::edges[e][0]], vertices[quad4::edges[e][1]], gid});

        if(corners.has_value()) {
          a2a[entity_kind<D>::corners].add_row(
            {(*corners)[quad4::edges[e][0]], (*corners)[quad4::edges[e][1]]});
        } // if
        a2a[interface_kind<D>()].add_row({interfaces[e]});
      } // for
      break;
    default:
      flog_fatal("invalid edge definition");
  } // switch

  return std::make_pair(std::move(entities), std::move(a2a));
} // create_cell_entities

// clang-format off
struct tet4 {
  static constexpr int num_edges{6};
  static constexpr std::uint32_t edges[num_edges][2] = {
    {0, 1}, /* v0, v1 */
    {1, 2},
    {2, 0},
    {0, 3},
    {1, 3},
    {2, 3}
  };

  static constexpr int num_faces{4};
  static constexpr std::uint32_t faces[num_faces][3] = {
    {0, 2, 1}, /* v0, v1, v2 */
    {0, 1, 3},
    {1, 2, 3},
    {2, 0, 3}
  };

  static constexpr int num_sides{12};
  static constexpr std::uint32_t sides[num_sides][3] = {
    {0, 2, 0}, /* v0, v1, f */
    {2, 1, 0},
    {1, 0, 0},

    {0, 1, 1},
    {1, 3, 1},
    {3, 0, 1},

    {1, 2, 2},
    {2, 3, 2},
    {3, 1, 2},

    {2, 0, 3},
    {0, 3, 3},
    {3, 2, 3}
  };
}; // struct tet4

struct hex8 {
  static constexpr int num_edges{12};
  static constexpr std::uint32_t edges[num_edges][2] = {
    {0, 1}, /* v0, v1 */
    {1, 2},
    {2, 3},
    {3, 0},
    {0, 4},
    {1, 5},
    {2, 6},
    {3, 7},
    {4, 5},
    {5, 6},
    {6, 7},
    {7, 4}
  };

  static constexpr int num_faces{6};
  static constexpr std::uint32_t faces[num_faces][4] = {
    {0, 3, 2, 1}, /* v0, v1, v2, v3 */
    {0, 1, 5, 4},
    {1, 2, 6, 5},
    {2, 3, 7, 6},
    {0, 4, 7, 3},
    {4, 5, 6, 7}
  };

  static constexpr int num_sides{24};
  static constexpr std::uint32_t sides[num_sides][3] = {
    {0, 3, 0}, /* v0, v1, f */
    {3, 2, 0},
    {2, 1, 0},
    {1, 0, 0},

    {0, 1, 1},
    {1, 5, 1},
    {5, 4, 1},
    {4, 0, 1},

    {1, 2, 2},
    {2, 6, 2},
    {6, 5, 2},
    {5, 1, 2},

    {2, 3, 3},
    {3, 7, 3},
    {7, 6, 3},
    {6, 2, 3},

    {0, 4, 4},
    {4, 7, 4},
    {7, 3, 4},
    {3, 0, 4},

    {4, 5, 5},
    {5, 6, 5},
    {6, 7, 5},
    {7, 4, 5}
  };
}; // struct hex8
// clang-format on

// 3D Edges.
template<std::size_t D, entity_kind<D> K>
auto
create_cell_entities(std::tuple<util::gid, util::id, util::id> const & cid,
  util::crs const & c2v,
  std::map<entity_kind<D>, util::crs> const &)
  -> std::enable_if_t<D == 3 && K == entity_kind<D>::edges,
    std::pair<util::crs, std::map<entity_kind<D>, util::crs>>> {
  auto [gid, id, cfaid] = cid;
  auto const & vertices = c2v[id];
  util::crs entities;
  std::map<entity_kind<D>, util::crs> a2a;

  switch(vertices.size()) {
    case 4 /* tet4 */:
      for(int e{0}; e < tet4::num_edges; ++e) {
        entities.add_row(
          {vertices[tet4::edges[e][0]], vertices[tet4::edges[e][1]]});
      } // for
      break;
    case 8 /* hex8 */:
      for(int e{0}; e < hex8::num_edges; ++e) {
        entities.add_row(
          {vertices[hex8::edges[e][0]], vertices[hex8::edges[e][1]]});
      } // for
      break;
    default:
      flog_fatal("invalid edge definition");
  } // switch

  return std::make_pair(std::move(entities), std::move(a2a));
} // create_cell_entities

// 3D Faces.
template<std::size_t D, entity_kind<D> K>
auto
create_cell_entities(std::tuple<util::gid, util::id, util::id> const & cid,
  util::crs const & c2v,
  std::map<entity_kind<D>, util::crs> const &)
  -> std::enable_if_t<D == 3 && K == entity_kind<D>::faces,
    std::pair<util::crs, std::map<entity_kind<D>, util::crs>>> {
  auto [gid, id, cfaid] = cid;
  auto const & vertices = c2v[id];
  util::crs entities;
  std::map<entity_kind<D>, util::crs> a2a;

  switch(vertices.size()) {
    case 4 /* tet4 */:
      for(int e{0}; e < tet4::num_faces; ++e) {
        entities.add_row({vertices[tet4::faces[e][0]],
          vertices[tet4::faces[e][1]],
          vertices[tet4::faces[e][2]]});
      } // for
      break;
    case 8 /* hex8 */:
      for(int e{0}; e < hex8::num_faces; ++e) {
        entities.add_row({vertices[hex8::faces[e][0]],
          vertices[hex8::faces[e][1]],
          vertices[hex8::faces[e][2]],
          vertices[hex8::faces[e][3]]});
      } // for
      break;
    default:
      flog_fatal("invalid face definition");
  } // switch

  return std::make_pair(std::move(entities), std::move(a2a));
} // create_cell_entities

// 3D Sides.
template<std::size_t D, entity_kind<D> K>
auto
create_cell_entities(std::tuple<util::gid, util::id, util::id> const & cid,
  util::crs const & c2v,
  std::map<entity_kind<D>, util::crs> const & aux)
  -> std::enable_if_t<D == 3 && K == entity_kind<D>::sides,
    std::pair<util::crs, std::map<entity_kind<D>, util::crs>>> {
  auto [gid, id, cfaid] = cid;
  auto const & vertices = c2v[id];
  util::crs entities;
  std::map<entity_kind<D>, util::crs> a2a;

  flog_assert(aux.count(entity_kind<D>::faces), "faces are required");

  auto const & faces = aux.at(entity_kind<D>::faces)[cfaid];

  auto const & corners =
    aux.count(entity_kind<D>::corners)
      ? std::make_optional(aux.at(entity_kind<D>::corners)[cfaid])
      : std::nullopt;

  // FIXME: Need to capture correct orientation information, i.e.,
  // invert vertices if ones complement set for face.
  switch(vertices.size()) {
    case 4 /* tet4 */:
      for(int s{0}; s < tet4::num_sides; ++s) {
        auto const fid = util::get_id(faces[tet4::sides[s][2]]);
        entities.add_row(
          {vertices[tet4::sides[s][0]], vertices[tet4::sides[s][1]], fid, gid});

        if(corners.has_value()) {
          a2a[entity_kind<D>::corners].add_row(
            {(*corners)[tet4::sides[s][0]], (*corners)[tet4::sides[s][1]]});
        } // if
        a2a[interface_kind<D>()].add_row({fid});
      } // for
      break;
    case 8 /* hex8 */:
      for(int s{0}; s < hex8::num_sides; ++s) {
        auto const fid = util::get_id(faces[hex8::sides[s][2]]);
        entities.add_row(
          {vertices[hex8::sides[s][0]], vertices[hex8::sides[s][1]], fid, gid});
        if(corners.has_value()) {
          a2a[entity_kind<D>::corners].add_row(
            {(*corners)[hex8::sides[s][0]], (*corners)[hex8::sides[s][1]]});
        } // if
        a2a[interface_kind<D>()].add_row({fid});
      } // for
      break;
    default:
      flog_fatal("invalid side definition");
  } // switch

  return std::make_pair(std::move(entities), std::move(a2a));
} // create_cell_entities

// 2D and 3D Corners.
template<std::size_t D, entity_kind<D> K>
auto
create_cell_entities(std::tuple<util::gid, util::id, util::id> const & cid,
  util::crs const & c2v,
  std::map<entity_kind<D>, util::crs> const &)
  -> std::enable_if_t<(D == 2 || D == 3) && K == entity_kind<D>::corners,
    std::pair<util::crs, std::map<entity_kind<D>, util::crs>>> {
  auto [gid, id, cfaid] = cid;
  auto const & vertices = c2v[id];
  util::crs entities;
  std::map<entity_kind<D>, util::crs> a2a;

  for(std::size_t v{0}; v < vertices.size(); ++v) {
    entities.add_row({vertices[v], gid});
  } // for

  return std::make_pair(std::move(entities), std::move(a2a));
} // create_cell_entities

} // namespace flsp::topo::unstructured::io

#endif // FLSP_TOPO_UNSTRUCTURED_MODELS_HH
