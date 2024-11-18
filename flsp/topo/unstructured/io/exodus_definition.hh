#ifndef FLSP_TOPO_UNSTRUCTURED_IO_EXODUS_DEFINITION_HH
#define FLSP_TOPO_UNSTRUCTURED_IO_EXODUS_DEFINITION_HH

/**
 * Exodus definition adapted from implementation by Marc Charest.
 */

#include "flsp/config.hh"
#include "flsp/topo/unstructured/io/definition_base.hh"
#include "flsp/topo/unstructured/io/types.hh"

#include <exodusII.h>

#include <iostream>
#include <memory>
#include <vector>

namespace flsp::topo::unstructured::io {

namespace detail {

template<class size, class block_t, auto is_fixed_block>
struct block_cursor {
  struct location {
    int block;
    size offset;
  };

  block_cursor(std::vector<std::pair<size, block_t>> && entities_per_block,
    size chunk_size)
    : base(0), block_counts(std::move(entities_per_block)),
      chunk_size(chunk_size) {}

  /**
   * check whether entity (id) is contained in current chunk.
   */
  bool contains(size id) const {
    return id >= base and id < next();
  }

  /**
   * Get id of the next entity to read.
   */
  size next() const {
    return (curr.size() == 0) ? 0 : base + curr.size();
  }

  void move(int block, size offset) {
    curr_loc = {block, offset};
    base = 0;
    for(int i = 0; i < block; i++) {
      base += block_counts[i].first;
    }
    base += offset;

    curr.clear();
  }

  location find_entity(size eid) {
    size curr = 0;
    location ret{0, 0};
    for(const auto & blkinfo : block_counts) {
      if(not is_fixed_block(blkinfo.second)) { // polyhedra
        curr += blkinfo.first;
        if(eid < curr)
          return ret;
        ++ret.block;
      }
      else { // fixed element so read by chunks
        if(eid < curr + blkinfo.first) { // first check if in block
          long long entities_remaining = blkinfo.first;
          ret.offset = 0;
          while(entities_remaining > 0) {
            entities_remaining -= chunk_size;
            if(eid < curr + ret.offset + chunk_size)
              return ret;
            ret.offset += chunk_size;
          }
        }
        curr += blkinfo.first;
        ++ret.block;
      }
    }

    // did not find entity (eid)
    ret.block = -1;
    ret.offset = -1;
    return ret;
  }

  void get(size id, std::vector<util::gid> & ret) {
    if(contains(id)) {
      auto sp = curr[id - base];
      ret.reserve(sp.size());
      std::copy(sp.begin(), sp.end(), std::back_inserter(ret));
    }
  }

  util::crs & data() {
    return curr;
  }

  const location & current_location() const {
    return curr_loc;
  }

  block_t get_block_type(size id) const {
    auto & binfo = block_counts.at(id);
    return binfo.second;
  }

  size get_base() const {
    return base;
  }

protected:
  size base;
  std::vector<std::pair<size, block_t>> block_counts;
  size chunk_size;
  util::crs curr;
  location curr_loc;
};

template<class size_t, class real>
class vertex_cursor
{
public:
  /**
   * \param num_verts global number of vertices to stream
   */
  vertex_cursor(size_t num_verts, size_t num_dims)
    : base_(0), num_verts_(num_verts), num_dims_(num_dims) {}

  /**
   * check whether entity (id) is contained in current chunk.
   */
  bool contains(size_t id) const {
    return id >= base_ and id < next();
  }

  /**
   * Get id of the next entity to read.
   */
  size_t next() const {
    return (size() == 0) ? 0 : base_ + size();
  }

  size_t size() const {
    return curr_.size() / num_dims_;
  }

  real * data() {
    return curr_.data();
  }

  const std::vector<real> & current() const {
    return curr_;
  }

  /**
   * get interval [begin, end) of entities in curr.
   */
  std::pair<size_t, size_t> interval() const {
    return std::make_pair(base_, base_ + size());
  }

  /**
   * Move cursor (prepare for inserting another chunk of entities).
   */
  void move_next(size_t chunksize) {
    base_ = next();
    if(base_ >= num_verts_)
      base_ = 0;
    curr_.clear();
    chunksize = std::min(chunksize, num_verts_ - base_);
    curr_.resize(chunksize * num_dims_);
  }

protected:
  size_t base_;
  size_t num_verts_;
  size_t num_dims_;
  std::vector<real> curr_;
};

/* tetra */
constexpr int tetra_sides = 4;
constexpr int tetra_size = 6;
constexpr int tetra_table[tetra_sides][tetra_size] = {
  /*      1              2               3               4 side */
  {1, 2, 4, 5, 9, 8},
  {2, 3, 4, 6, 10, 9},
  {1, 4, 3, 8, 10, 7},
  {1, 3, 2, 7, 6, 5} /* nodes  */
};

/* hex */
constexpr int hex_sides = 6;
constexpr int hex_size = 9;
constexpr int hex_table[hex_sides][hex_size] = {
  /*         1                        2 side   */
  {1, 2, 6, 5, 9, 14, 17, 13, 26},
  {2, 3, 7, 6, 10, 15, 18, 14, 25}, /* nodes  */
  /*         3                        4 side   */
  {3, 4, 8, 7, 11, 16, 19, 15, 27},
  {1, 5, 8, 4, 13, 20, 16, 12, 24}, /* nodes  */
  /*         5                        6 side   */
  {1, 4, 3, 2, 12, 11, 10, 9, 22},
  {5, 6, 7, 8, 17, 18, 19, 20, 23} /* nodes  */
};

//! An enumeration to keep track of element types
enum class block_t {
  tri,
  quad,
  polygon,
  tet,
  hex,
  polyhedron,
  unknown,
  empty,
  invalid
};

template<typename U>
struct block_stats_t {
  U num_elem_this_blk = 0;
  U num_faces_per_elem = 0;
  U num_edges_per_elem = 0;
  U num_nodes_per_elem = 0;
  U num_attr = 0;
  char elem_type[MAX_STR_LENGTH];
};

inline bool
is_fixed_block(block_t blk) {
  return not(blk == block_t::polygon or blk == block_t::polyhedron);
}

using bcursor = block_cursor<std::size_t, block_t, is_fixed_block>;
using vcursor = vertex_cursor<std::size_t, double>;

template<class U>
block_t
get_block_type(const block_stats_t<U> & stats) {
  if(strcasecmp("nsided", stats.elem_type) == 0)
    return block_t::polygon;
  else if(strcasecmp("nfaced", stats.elem_type) == 0)
    return block_t::polyhedron;
  else if(strcasecmp("tri", stats.elem_type) == 0 ||
          strcasecmp("tri3", stats.elem_type) == 0)
    return block_t::tri;
  else if(strcasecmp("quad", stats.elem_type) == 0 ||
          strcasecmp("quad4", stats.elem_type) == 0 ||
          strcasecmp("shell", stats.elem_type) == 0 ||
          strcasecmp("shell4", stats.elem_type) == 0)
    return block_t::quad;
  else if(strcasecmp("tet", stats.elem_type) == 0 ||
          strcasecmp("tetra", stats.elem_type) == 0 ||
          strcasecmp("tet4", stats.elem_type) == 0)
    return block_t::tet;
  else if(strcasecmp("hex", stats.elem_type) == 0 ||
          strcasecmp("hex8", stats.elem_type) == 0)
    return block_t::hex;
  else
    return block_t::unknown;
}

int
open(const std::string & name, std::ios_base::openmode mode) {
#ifdef DEBUG
  ex_opts(EX_ABORT | EX_VERBOSE);
#endif

  // size of floating point variables used in app.
  int app_word_size = sizeof(double);

  if((mode & std::ios_base::in) == std::ios_base::in) {
    // size of floating point stored in name.
    int exo_word_size = 0;
    // the version number
    float version;

    // open the file
    auto exoid =
      ex_open(name.c_str(), EX_READ, &app_word_size, &exo_word_size, &version);
    if(exoid < 0)
      flog_fatal("Problem opening exodus file, ex_open() returned " << exoid);

    // This sets the file to read IDs as 64 bit.  If the file does not have
    // 64 bit IDs, it should have no effect.
    ex_set_int64_status(exoid, EX_ALL_INT64_API);

    return exoid;
  }
  else if((mode & std::ios_base::out) == std::ios_base::out) {
    // size of floating point to be stored in file.
    // change to float to save space
    int exo_word_size = sizeof(double);

    // determine the file creation mode
    int cmode = (mode & std::ios_base::app) == std::ios_base::app ? EX_NOCLOBBER
                                                                  : EX_CLOBBER;

    // create file
    auto exoid = ex_create(name.c_str(), cmode, &app_word_size, &exo_word_size);
    if(exoid < 0)
      flog_fatal("Problem writing exodus file, ex_create() returned " << exoid);

    return exoid;
  }
  else {

    flog_fatal("Unknown file mode");
    return -1;
  }
}

inline void
close(int exoid) {
  auto status = ex_close(exoid);
  if(status)
    flog_fatal("Problem closing exodus file, ex_close() returned " << exoid);
}

inline ex_init_params
read_params(int exoid, std::size_t dimension) {
  ex_init_params exo_params;
  auto status = ex_get_init_ext(exoid, &exo_params);
  if(status) {
    flog_fatal(
      "Problem getting exodus file parameters, ex_get_init_ext() returned "
      << status);
  }

  if(dimension != std::size_t(exo_params.num_dim)) {
    flog_fatal("Exodus dimension mismatch: Expected dimension ("
               << dimension << ") != Exodus dimension (" << exo_params.num_dim
               << ")");
  }

  return exo_params;
}

inline void
write_params(int exoid,
  std::size_t dimension,
  const ex_init_params & exo_params) {
  if(dimension != std::size_t(exo_params.num_dim))
    flog_fatal("Exodus dimension mismatch: Expected dimension ("
               << dimension << ") /= Exodus dimension (" << exo_params.num_dim
               << ")");

  auto status = ex_put_init_ext(exoid, &exo_params);
  if(status)
    flog_fatal(
      "Problem putting exodus file parameters, ex_put_init_ext() returned "
      << status);
}

inline bool
is_int64(int exoid) {
  return (ex_int64_status(exoid) & EX_IDS_INT64_API);
}

template<typename U>
void
read_block_stats(int exoid_,
  ex_entity_id blk_id,
  ex_entity_type entity_type,
  block_stats_t<U> & block_stats) {
  auto status = ex_get_block(exoid_,
    entity_type,
    blk_id,
    block_stats.elem_type,
    &block_stats.num_elem_this_blk,
    &block_stats.num_nodes_per_elem,
    &block_stats.num_edges_per_elem,
    &block_stats.num_faces_per_elem,
    &block_stats.num_attr);
  if(status)
    flog_fatal("Problem reading block, ex_get_block() returned " << status);
}

template<class U>
std::vector<std::size_t>
read_block_ids(int exoid, ex_entity_type obj_type, std::size_t num_blocks) {
  using ex_index_t = U;

  std::vector<std::size_t> ids(num_blocks);

  if(num_blocks > 0) {
    std::vector<ex_index_t> block_ids(num_blocks);
    auto status = ex_get_ids(exoid, obj_type, block_ids.data());
    if(status)
      flog_fatal("Problem reading block ids, ex_get_ids() returned " << status);

    // now convert them
    std::transform(
      block_ids.begin(), block_ids.end(), ids.begin(), [](auto id) {
        return id;
      });
  }

  return ids;
}

static constexpr std::size_t CHUNK_SIZE = 256;

template<class U>
block_t
read_block(int exoid,
  ex_entity_id blkid,
  ex_entity_type entity_type,
  bcursor & blk_cursor,
  U start) {
  using stats_t = block_stats_t<U>;
  stats_t stats;
  read_block_stats(exoid, blkid, EX_ELEM_BLOCK, stats);

  if(not stats.num_elem_this_blk)
    return block_t::empty;

  std::vector<int> counts;
  std::vector<U> indices;

  block_t ret;

  // polygon data
  if(strcasecmp("nsided", stats.elem_type) == 0) {
    // the number of nodes per element is really the number of nodes
    // in the whole block
    auto num_nodes_this_blk = stats.num_nodes_per_elem;

    // get the number of nodes per element
    counts.resize(stats.num_elem_this_blk);
    auto status = ex_get_entity_count_per_polyhedra(
      exoid, entity_type, blkid, counts.data());
    if(status)
      flog_fatal("Problem getting element node numbers, "
                 << "ex_get_entity_count_per_polyhedra() returned " << status);

    // read element definitions
    indices.resize(num_nodes_this_blk);
    status =
      ex_get_conn(exoid, entity_type, blkid, indices.data(), nullptr, nullptr);
    if(status)
      flog_fatal("Problem getting element connectivity, ex_get_elem_conn() "
                 << "returned " << status);

    ret = block_t::polygon;
  }
  // polygon data
  else if(strcasecmp("nfaced", stats.elem_type) == 0) {
    // the number of faces per element is really the number of
    // faces in the whole block ( includes duplicate / overlapping
    // nodes )
    auto num_face_this_blk = stats.num_faces_per_elem;

    // get the number of nodes per element
    counts.resize(stats.num_elem_this_blk);
    auto status = ex_get_entity_count_per_polyhedra(
      exoid, entity_type, blkid, counts.data());
    if(status)
      flog_fatal("Problem reading element node info, "
                 << "ex_get_entity_count_per_polyhedra() returned " << status);

    // read element definitions
    indices.resize(num_face_this_blk);
    status =
      ex_get_conn(exoid, entity_type, blkid, nullptr, nullptr, indices.data());
    if(status)
      flog_fatal("Problem getting element connectivity, ex_get_conn() "
                 << "returned " << status);

    ret = block_t::polyhedron;
  }
  // fixed element size
  else {
    U end = start + CHUNK_SIZE;
    end = std::min(end, stats.num_elem_this_blk);
    auto num_elem_this_blk = end - start;

    // set the counts
    counts.resize(num_elem_this_blk);
    std::fill(counts.begin(), counts.end(), stats.num_nodes_per_elem);

    // read element definitions
    indices.resize(num_elem_this_blk * stats.num_nodes_per_elem);
    auto status = ex_get_partial_conn(exoid,
      EX_ELEM_BLOCK,
      blkid,
      start + 1,
      end - start,
      indices.data(),
      0,
      0);
    if(status)
      flog_fatal("Problem getting element connectivity, ex_get_elem_conn() "
                 << "returned " << status);

    // return element type
    if(strcasecmp("tri", stats.elem_type) == 0 ||
       strcasecmp("tri3", stats.elem_type) == 0)
      ret = block_t::tri;
    else if(strcasecmp("quad", stats.elem_type) == 0 ||
            strcasecmp("quad4", stats.elem_type) == 0 ||
            strcasecmp("shell", stats.elem_type) == 0 ||
            strcasecmp("shell4", stats.elem_type) == 0)
      ret = block_t::quad;
    else if(strcasecmp("tet", stats.elem_type) == 0 ||
            strcasecmp("tetra", stats.elem_type) == 0 ||
            strcasecmp("tet4", stats.elem_type) == 0)
      ret = block_t::tet;
    else if(strcasecmp("hex", stats.elem_type) == 0 ||
            strcasecmp("hex8", stats.elem_type) == 0)
      ret = block_t::hex;
    else {
      flog_fatal("Unknown block type, " << stats.elem_type);
      ret = block_t::unknown;
    }
  } // element type

  { // filter block

    // create cells in mesh
    std::size_t base = 0;
    for(std::size_t e = 0; e < counts.size(); e++) {
      // get the number of nodes
      auto cnt = counts[e];
      // make vertex indices zero-based and copy local vertices into a new row
      blk_cursor.data().add_row(flecsi::util::transform_view(
        flecsi::util::span(indices.data() + base, cnt),
        [](auto x) { return x - 1; }));
      base += cnt;
    }
  }

  return ret;
}

} // namespace detail
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

#if defined(FLECSI_SP_ENABLE_EXODUSII)
const inline bool register_exodus_1d_ =
  io_factory<1>::instance().register_type("exo", exodus_handler<1>);
const inline bool register_exodus_2d_ =
  io_factory<2>::instance().register_type("exo", exodus_handler<2>);
const inline bool register_exodus_3d_ =
  io_factory<3>::instance().register_type("exo", exodus_handler<3>);
#endif

} // namespace flsp::topo::unstructured::io

#endif // FLSP_TOPO_UNSTRUCTURED_IO_EXODUS_DEFINITION_HH
