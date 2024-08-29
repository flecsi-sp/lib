#ifndef FLSP_TOPO_UNSTRUCTURED_IO_EXODUS_DEFINITION_HH
#define FLSP_TOPO_UNSTRUCTURED_IO_EXODUS_DEFINITION_HH

/**
 * Exodus definition adapted from implementation by Marc Charest.
 */

#include "flsp/topo/unstructured/io/definition_base.hh"
#include "flsp/topo/unstructured/io/detail/exodus.hh"
#include "flsp/topo/unstructured/io/types.hh"

#include <exodusII.h>

#include <iostream>
#include <memory>
#include <vector>

namespace flsp::topo::unstructured::io {

/*----------------------------------------------------------------------------*
 * Base
 *----------------------------------------------------------------------------*/

template<std::size_t D>
struct exodus_definition : definition_base<D> {
  using index = std::size_t;
  static constexpr std::size_t CHUNK_SIZE = 256;
  using id = std::size_t;

  exodus_definition(const std::string & filename) {
    exoid_ = detail::open(filename, std::ios_base::in);
    if(exoid_ < 0)
      flog_fatal("Problem reading exodus file");

    exo_params = detail::read_params(exoid_, D);
    vert_cursor = std::make_unique<detail::vcursor>(exo_params.num_nodes, D);

    auto num_elem_blk = exo_params.num_elem_blk;

    if(detail::is_int64(exoid_))
      elem_blk_ids =
        detail::read_block_ids<long long>(exoid_, EX_ELEM_BLOCK, num_elem_blk);
    else
      elem_blk_ids =
        detail::read_block_ids<int>(exoid_, EX_ELEM_BLOCK, num_elem_blk);

    elem_blk_types.resize(elem_blk_ids.size());

    // read block counts to initialize block_cursor (so it can find what block
    // contains an entity)
    std::vector<std::pair<std::size_t, detail::block_t>> block_counts;
    for(auto blkid : elem_blk_ids) {
      if(detail::is_int64(exoid_)) {
        detail::block_stats_t<long long> stats;
        read_block_stats(exoid_, blkid, EX_ELEM_BLOCK, stats);
        block_counts.emplace_back(
          stats.num_elem_this_blk, get_block_type(stats));
      }
      else {
        detail::block_stats_t<int> stats;
        read_block_stats(exoid_, blkid, EX_ELEM_BLOCK, stats);
        block_counts.emplace_back(
          stats.num_elem_this_blk, get_block_type(stats));
      }
    }
    blk_cursor =
      std::make_unique<detail::bcursor>(std::move(block_counts), CHUNK_SIZE);
  }

  ~exodus_definition() {
    close(this->exoid_);
  }

  exodus_definition(const exodus_definition &) = delete;
  exodus_definition & operator=(const exodus_definition &) = delete;

  std::tuple<util::crs, std::optional<face_info>, std::optional<mat_ids>>
  cell_data(iota_view const & r) const override {
    util::crs p2v;
    for(auto p : r) {
      std::vector<util::gid> def;
      stream(p, def);
      p2v.add_row(def);
    }
    return std::make_tuple(p2v, std::nullopt, std::nullopt);
  }

  template<class VID>
  void build_intermediary_from_vertices(flecsi::Dimension idim,
    std::size_t,
    const std::vector<VID> & verts,
    util::crs & inter) const;

  util::gid num_entities(entity_kind<D> k) const override {
    flog_assert(k == entity_kind<D>::cells || k == entity_kind<D>::vertices,
      "invalid entity_kind");
    return k == entity_kind<D>::vertices ? exo_params.num_nodes
                                         : exo_params.num_elem;
  }

  std::tuple<std::vector<util::point<D>>, std::optional<bnd_ids>> vertex_data(
    iota_view const & r) const override {
    std::vector<util::point<D>> c;
    for(auto v : r) {
      while(not vert_cursor->contains(v)) {
        read_next_vertices();
      }
      auto [beg, end] = vert_cursor->interval();
      auto num_nodes = end - beg;
      const auto & vertices = vert_cursor->current();
      util::point<D> vert;
      for(std::size_t i = 0; i < D; i++) {
        vert[i] = vertices[i * num_nodes + (v - beg)];
      }
      c.push_back(vert);
    } // for

    return std::tuple(std::move(c), std::nullopt);
  }

protected:
  detail::block_t stream_block(int block, std::size_t offset) const {
    if(std::size_t(block) >= elem_blk_ids.size())
      return detail::block_t::invalid;
    detail::block_t blktype;

    auto blkid = elem_blk_ids[block];

    blk_cursor->move(block, offset);
    if(detail::is_int64(exoid_))
      blktype = detail::read_block<long long>(
        exoid_, blkid, EX_ELEM_BLOCK, *blk_cursor, offset);
    else
      blktype = detail::read_block<int>(
        exoid_, blkid, EX_ELEM_BLOCK, *blk_cursor, offset);

    elem_blk_types[block] = blktype;

    return blktype;
  }

  void stream(index entity_id, std::vector<util::gid> & ret) const {
    if(not stream_contains(entity_id)) {
      auto loc = blk_cursor->find_entity(entity_id);
      if(loc.block == -1) {
        flog_fatal("Problem finding entity in file: " << entity_id);
      }
      stream_block(loc.block, loc.offset);
    }

    blk_cursor->get(entity_id, ret);
  }

  bool stream_contains(index entity_id) const {
    return blk_cursor->contains(entity_id);
  }

  std::vector<double> read_point_coords(std::size_t num_nodes) {
    if(num_nodes <= 0)
      flog_fatal(
        "Exodus file has zero nodes, or parmeters haven't been read yet.");

    // read nodes
    std::vector<double> vertex_coord(D * num_nodes);

    // exodus is kind enough to fetch the data in the double type we ask for
    auto status = ex_get_coord(exoid_,
      vertex_coord.data(),
      vertex_coord.data() + num_nodes,
      vertex_coord.data() + 2 * num_nodes);

    if(status)
      flog_fatal("Problem getting vertex coordinates from exodus file, "
                 << " ex_get_coord() returned " << status);

    return vertex_coord;
  }

  void write_point_coords(int exo_id,
    const std::vector<double> & vertex_coord) {
    if(vertex_coord.empty())
      return;

    auto num_nodes = vertex_coord.size() / D;

    // exodus is kind enough to fetch the data in the double type we ask for
    auto status = ex_put_coord(exo_id,
      vertex_coord.data(),
      vertex_coord.data() + num_nodes,
      vertex_coord.data() + 2 * num_nodes);

    if(status)
      flog_fatal("Problem putting vertex coordinates to exodus file, "
                 << " ex_put_coord() returned " << status);
  }

  ex_init_params get_params() const {
    return exo_params;
  }

  void read_next_vertices() const {
    vert_cursor->move_next(CHUNK_SIZE);
    auto [beg, end] = vert_cursor->interval();
    auto num_nodes = end - beg;
    auto status = ex_get_partial_coord(exoid_,
      beg + 1,
      num_nodes,
      vert_cursor->data(),
      vert_cursor->data() + num_nodes,
      vert_cursor->data() + 2 * num_nodes);

    if(status)
      flog_fatal("Problem getting vertex coordinates from exodus file, "
                 << " ex_get_partial_coord() returned " << status);
  }

  const detail::bcursor & get_block_cursor() const {
    return *blk_cursor;
  }

  index block_id(std::size_t ind) const {
    return elem_blk_ids[ind];
  }

  detail::block_t block_type(std::size_t ind) const {
    return elem_blk_types[ind];
  }

protected:
  int exoid_;
  ex_init_params exo_params;
  std::vector<index> elem_blk_ids;
  mutable std::vector<detail::block_t> elem_blk_types;
  bool int64;
  mutable std::unique_ptr<detail::bcursor> blk_cursor;
  mutable std::unique_ptr<detail::vcursor> vert_cursor;
}; // struct exodus_definition

/*----------------------------------------------------------------------------*
 * 1D
 *----------------------------------------------------------------------------*/

template<>
template<class VID>
void
exodus_definition<1>::build_intermediary_from_vertices(flecsi::Dimension idim,
  std::size_t,
  const std::vector<VID> & verts,
  util::crs & inter) const {
  flog_assert(idim == 1, "Invalid dimension: " << idim);
  for(auto v0 = verts.begin(), v1 = std::next(v0); v0 != verts.end();
      ++v0, ++v1) {
    if(v1 == verts.end())
      v1 = verts.begin();
    inter.add_row({*v0, *v1});
  }
}

/*----------------------------------------------------------------------------*
 * 2D
 *----------------------------------------------------------------------------*/

template<>
template<class VID>
void
exodus_definition<2>::build_intermediary_from_vertices(flecsi::Dimension idim,
  std::size_t,
  const std::vector<VID> & verts,
  util::crs & inter) const {
  if(idim == 1) { // edges and faces
    for(auto v0 = verts.begin(), v1 = std::next(v0); v0 != verts.end();
        ++v0, ++v1) {
      if(v1 == verts.end())
        v1 = verts.begin();
      inter.add_row({*v0, *v1});
    }
  }
  else if(idim == 4 /*D + 2*/) { // corners
    for(auto v : verts) {
      inter.add_row({v});
    }
  }
  else {
    flog_fatal("Invalid dimension: " << idim);
  }
}

/*----------------------------------------------------------------------------*
 * 3D
 *----------------------------------------------------------------------------*/

template<>
template<class VID>
void
exodus_definition<3>::build_intermediary_from_vertices(flecsi::Dimension idim,
  std::size_t cell,
  const std::vector<VID> & verts,
  util::crs & inter) const {
  if(idim == 1 || idim == 2) { // edges or faces
    // TODO: make more efficient (reduce amount of block lookups)
    auto loc = blk_cursor->find_entity(cell);
    auto btype = blk_cursor->get_block_type(loc.block);
    const int * side_indices = nullptr;
    int num_sides = 0;
    int num_side_points = 0;
    int side_size = 0;
    if(btype == detail::block_t::hex) {
      side_indices = &detail::hex_table[0][0];
      num_sides = detail::hex_sides;
      side_size = detail::hex_size;
      num_side_points = 4;
    }
    else if(btype == detail::block_t::tet) {
      side_indices = &detail::tetra_table[0][0];
      num_sides = detail::tetra_sides;
      side_size = detail::tetra_size;
      num_side_points = 3;
    }
    std::vector<index> temp_vs(num_side_points);

    for(int i = 0; i < num_sides; ++i) {
      for(int j = 0; j < num_side_points; ++j) {
        auto id = side_indices[i * side_size + j] - 1;
        temp_vs[j] = verts[id];
      }
      if(idim == 2) {
        // faces
        inter.add_row(temp_vs.begin(), temp_vs.end());
      }
      else {
        // edges
        for(auto v0 = std::prev(temp_vs.end()), v1 = temp_vs.begin();
            v1 != temp_vs.end();
            v0 = v1, ++v1)
          inter.add_row({*v0, *v1});
      }
    }
  }
  else if(idim == 5 /*D + 2*/) { // corners
    // One corner for each vertex of the cell
    for(auto v : verts) {
      inter.add_row({v});
    }
  }
  else {
    flog_fatal("Invalid dimension: " << idim);
  }
}

template<std::size_t D>
std::unique_ptr<definition_base<D>>
exodus_handler(const std::string & fname,
  std::optional<std::vector<std::string>>,
  std::optional<std::vector<std::string>>,
  MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  if(rank == 0) {
    return std::make_unique<exodus_definition<D>>(fname);
  }
  else {
    return std::make_unique<undefined_definition<D>>();
  } // if
}
const inline bool register_exodus_1d_ =
  io_factory<1>::instance().register_type("exo", exodus_handler<1>);
const inline bool register_exodus_2d_ =
  io_factory<2>::instance().register_type("exo", exodus_handler<2>);
const inline bool register_exodus_3d_ =
  io_factory<3>::instance().register_type("exo", exodus_handler<3>);

} // namespace flsp::topo::unstructured::io

#endif // FLSP_TOPO_UNSTRUCTURED_IO_EXODUS_DEFINITION_HH
