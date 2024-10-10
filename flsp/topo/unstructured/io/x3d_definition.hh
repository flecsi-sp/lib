#ifndef FLSP_TOPO_UNSTRUCTURED_IO_X3D_DEFINITION_H
#define FLSP_TOPO_UNSTRUCTURED_IO_X3D_DEFINITION_H

#include "flsp/topo/unstructured/io/definition_base.hh"
#include "flsp/topo/unstructured/io/types.hh"
#include "flsp/topo/unstructured/util/common.hh"

#include <flecsi/flog.hh>
#include <flecsi/util/common.hh>

#include <cmath>
#include <cstddef>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>

namespace flsp::topo::unstructured::io {

static void
assert_str(std::string_view expect, std::string_view got) {
  flog_assert(expect == got,
    "Error parsing X3D file: expected '" << expect << "' but received '" << got
                                         << ";");
}

struct x3d_face {
  std::size_t id;
  std::size_t num_nodes;
  std::vector<std::size_t> nodes;
  std::size_t neighbor;
};

struct x3d_cell {
  std::size_t id;
  std::size_t num_faces;
  std::vector<std::size_t> faces;
};

struct x3d_header {
  int process; // processor for which this file is targetted (1-based)
  std::size_t numdim; // spatial dimension
  int materials; // total number of materials in the mesh
  // the remaining entries are all for the current processor
  std::size_t nodes; // number of logically distinct nodes
  std::size_t faces; // number of logically distinct faces
  std::size_t elements; // number of cells
  std::size_t ghost_nodes; // number of ghost nodes
  std::size_t slaved_nodes; // number of slaved nodes
  // number of master nodes for most constrained slave
  // this will typically be 2 or 4 (depending on numdim), but should not be
  // zero even if there are no slaved nodes.
  int nodes_per_slave;
  int nodes_per_face; // number of nodes needed for the most complex face
  int faces_per_cell; // number of nodes needed for the most complex cell
  int node_data_fields; // number of node-centered data fields for simulation
  int cell_data_fields; // number of cell-centered data fields for simulation

  static constexpr auto X3D_TOKEN = "x3dtoflag ascii";

  void read(std::fstream & fh) {
    fh.seekg(0);
    std::string tok;

    std::getline(fh, tok);
    assert_str(tok, X3D_TOKEN);

    std::getline(fh, tok);
    assert_str(tok, "header");

    fh >> tok >> process;
    fh >> tok >> numdim;
    fh >> tok >> materials;
    fh >> tok >> nodes;
    fh >> tok >> faces;
    fh >> tok >> elements;
    fh >> tok >> ghost_nodes;
    fh >> tok >> slaved_nodes;
    fh >> tok >> nodes_per_slave;
    fh >> tok >> nodes_per_face;
    fh >> tok >> faces_per_cell;
    fh >> tok >> node_data_fields;
    fh >> tok >> cell_data_fields;

    std::getline(fh, tok); // newline from previous line
    std::getline(fh, tok);
    assert_str(tok, "end_header");
  } // read
};

struct x3d_seek {
  static constexpr std::size_t node_width = 80;

  x3d_seek(std::fstream & fh, const x3d_header & h) : cell_pos{0}, header{h} {
    seek_nodes(fh);
  }

  void seek_nodes(std::fstream & fh) {
    fh.seekg(0);
    std::string tok;
    do {
      std::getline(fh, tok);
    } while(tok != "nodes");
    node_pos = fh.tellg();
  }

  constexpr std::size_t face_pos() const {
    return node_pos + node_width * header.nodes + 10; // end_nodes
  }

  std::size_t node_pos, cell_pos;
  const x3d_header & header;
};

// Utility functions for the boundary and region readers

// The logic of the streaming approach requires that you attempt to read the
// value _after_ the last value and return something:
// - Loop over all elements in a given range.
// - If the current element in the reader has the same index as the current
//   element from the range, process that element from the reader and then
//   advance the reader to load the next element.
// - Once we process the last element from the reader, you "load" the element
//   _after_ the last element.
// - We never actually use that one-past-the-last element, because it will have
//   an element index of std::numeric_limits::max, which never matches the
//   index of the current element in the range.
// Therefore you cannot error out if you read a "bad" element, because that's
// part of normal operation.  Thus bad_value and bad_element exist to create
// fake "elements" after the last real element, and the fake elements have to be
// set up so that they will never be parsed as real elements.
template<typename T>
T
bad_value() {
  return std::numeric_limits<T>::max();
}

template<typename T>
T
read_next_value(std::ifstream & fin) {
  std::string line;
  std::getline(fin, line);
  if(line.empty()) {
    // See explanatory note for bad_value
    return bad_value<T>();
  }
  else {
    std::stringstream ss(line);
    T value;
    ss >> value;
    return value;
  }
}

// Support class for reading boundary files (a.k.a. point subset files)

class boundary_file_reader
{
private:
  using index_t = std::size_t;

  std::ifstream id_stream;

  auto next_element() {
    // change to 0-indexed
    return read_next_value<index_t>(id_stream) - 1;
  }

  auto bad_element() {
    // See explanatory note for bad_value
    return bad_value<index_t>();
  }

public:
  using element_t = index_t;

  element_t head;

  boundary_file_reader(std::string const & filename) : id_stream(filename) {
    // Load first element
    head = next_element();
    if(!*this) {
      flog_fatal("empty boundary file \"" << filename << "\"");
    }
  }

  void advance() {
    head = next_element();
    if(!*this) {
      head = bad_element();
    }
  }

  explicit operator bool() const {
    return bool(id_stream);
  }
};

// Support class for reading region files (a.k.a. cell subset files)

class region_file_reader
{
private:
  using index_t = std::size_t;

  std::ifstream id_stream;
  std::ifstream vf_stream;
  bool has_volume_fractions;

  auto next_element() {
    return std::make_pair(
      read_next_value<index_t>(id_stream) - 1, // change to 0-indexed
      has_volume_fractions ? read_next_value<double>(vf_stream) : 1.0);
  }

  auto bad_element() {
    // See explanatory note for bad_value
    return std::make_pair(bad_value<index_t>(), bad_value<double>());
  }

public:
  using element_t = std::pair<index_t, double>;

  element_t head;

  region_file_reader(std::string const & filename)
    : id_stream(filename), vf_stream(filename), has_volume_fractions{false} {
    // Advance vf_stream past IDs, count IDs
    std::string line;
    index_t id_count{0};
    while(std::getline(vf_stream, line)) {
      if(line.empty()) {
        has_volume_fractions = true;
        break;
      }
      ++id_count;
    }
    // Count volume fractions and verify matching counts
    if(has_volume_fractions) {
      auto vf_pos = vf_stream.tellg();
      index_t vf_count{0};
      while(std::getline(vf_stream, line)) {
        ++vf_count;
      }
      assert(vf_count == id_count);
      vf_stream.clear();
      vf_stream.seekg(vf_pos, std::ios::beg);
    }
    // Load first element
    head = next_element();
    if(!*this) {
      flog_fatal("empty region file \"" << filename << "\"");
    }
  }

  void advance() {
    head = next_element();
    if(!*this) {
      head = bad_element();
    }
  }

  explicit operator bool() const {
    return has_volume_fractions ? bool(vf_stream) : bool(id_stream);
  }
};

template<std::size_t D>
struct x3d_definition : definition_base<D> {

  static_assert(required_keys<entity_kind<D>,
                  entity_kind<D>::cells,
                  entity_kind<D>::vertices,
                  entity_kind<D>::faces>::value,
    "required entity kind undefined");

  using size = std::size_t;

  static constexpr std::size_t num_dims = D;
  static constexpr flecsi::Dimension dimension() {
    return D;
  }

  x3d_definition(const std::string & fname,
    const std::optional<std::vector<std::string>> & matfiles,
    const std::optional<std::vector<std::string>> & bndfiles)
    : fh(fname, std::ios_base::in) {
    flog_assert(fh.is_open(), "Error opening file: " << fname);

    header.read(fh);

    flog_assert(
      header.process == 1, "X3D reader does not support distributed X3D files");

    flog_assert(header.numdim == D,
      "X3D: expected " << D << " dimensions, but got " << header.numdim);

    seeker.emplace(fh, header);

    if(matfiles.has_value()) {
      for(const std::string & matfname : *matfiles) {
        matreaders.emplace_back(matfname.c_str());
      }
    }

    if(bndfiles.has_value()) {
      for(const std::string & bndfname : *bndfiles) {
        bndreaders.emplace_back(bndfname.c_str());
      }
    }

    read_all_nodes();
    read_all_faces();
    read_all_cells();
  }

  virtual std::
    tuple<util::crs, std::optional<face_info>, std::optional<mat_ids>>
    cell_data(const iota_view & rng) const override {
    std::vector<size> faces;
    face_info finfo;
    auto & c2f_proc = finfo.c2f;
    auto & f2v_proc = finfo.f2v;
    auto & p2m = finfo.p2m; // p2m for f2v (faces)
    util::crs c2v_proc;

    c2f_proc.offsets.reserve(rng.size() + 1);
    c2f_proc.values.reserve(
      rng.size() * c2f[0].size()); // could have better estimate
    faces.reserve(
      0.5 * rng.size() * c2f[0].size()); // could have better estimate
    c2v_proc.offsets.reserve(rng.size() + 1);

    std::optional<mat_ids> materials;
    // first check if we will be reading any material data (to avoid logic if
    // not)
    for(const auto & mread : matreaders) {
      // assuming material files and ranges are monotonically increasing
      if(mread.head.first <= rng.back()) {
        if(!materials.has_value()) {
          materials.emplace();
          (*materials).resize(matreaders.size());
        }
      }
    }
    for(auto i : rng) {
      if(materials.has_value())
        check_materials(*materials, i);

      std::vector<size> cell_faces;
      std::vector<size> cell_verts;
      auto c2f_row = c2f[i];
      // temporary assume fixed elements based on number of faces (to set c2v
      // according to element model).
      bool is_hex = (D == 3 && c2f_row.size() == 6);
      bool is_quad_or_hex =
        ((D == 2 && c2f_row.size() == 4) || is_hex) ? true : false;
      for(auto f : c2f_row) {
        size face = util::get_id(f);

        faces.push_back(face);
        cell_faces.push_back(f);

        if(is_quad_or_hex) { // use element model for node ordering
          std::vector<size> face_verts;
          bool all_new = true;
          for(size vi : f2v[face]) {
            auto it = std::find(cell_verts.begin(), cell_verts.end(), vi);
            if(it == cell_verts.end())
              face_verts.push_back(vi);
            else
              all_new = false;
          }
          if(all_new) {
            // order vertices based on face orientation
            if(util::sign_bit(f))
              std::reverse(face_verts.begin(), face_verts.end());
            if(is_hex && cell_verts.size() > 0) {
              /* if hex then reverse second face orientation and align with
                 previously added face
              */
              std::reverse(face_verts.begin(), face_verts.end());
              std::size_t starting_pos = 0;
              auto face_count = [&](size v1, size v2) {
                int cnt = 0;
                for(auto ff : c2f_row) {
                  size ffid = util::get_id(ff);
                  auto it1 = std::find(f2v[ffid].begin(), f2v[ffid].end(), v1);
                  auto it2 = std::find(f2v[ffid].begin(), f2v[ffid].end(), v2);
                  if(it1 != f2v[ffid].end() && it2 != f2v[ffid].end())
                    ++cnt;
                }
                return cnt;
              };
              // alignment node shares 2 faces with vertex on previously added
              // face
              for(std::size_t ii = 0; ii < face_verts.size(); ++ii) {
                if(face_count(cell_verts[0], face_verts[ii]) == 2)
                  starting_pos = ii;
              }
              if(starting_pos > 0) {
                std::rotate(face_verts.begin(),
                  face_verts.begin() + starting_pos,
                  face_verts.end());
              }
            }
            cell_verts.insert(
              cell_verts.end(), face_verts.begin(), face_verts.end());
          }
        }
        else {
          for(size vi : f2v[face]) {
            auto it = std::find(cell_verts.begin(), cell_verts.end(), vi);
            if(it == cell_verts.end())
              cell_verts.push_back(vi);
          }
        }
      }
      c2f_proc.add_row(cell_faces);
      c2v_proc.add_row(cell_verts);
    }

    flecsi::util::force_unique(faces);

    p2m.reserve(faces.size());
    f2v_proc.offsets.reserve(faces.size() + 1);
    f2v_proc.values.reserve(
      faces.size() * f2v[0].size()); // could have better estimate

    for(auto f : faces) {
      std::vector<size> verts;
      p2m.push_back(f);
      for(auto v : f2v[f])
        verts.push_back(v);
      f2v_proc.add_row(verts);
    }

    return std::tuple(
      std::move(c2v_proc), std::move(finfo), std::move(materials));
  }

protected:
  void check_materials(mat_ids & matids, util::gid i) const {
    flog_assert(matids.size() == matreaders.size(),
      "Incorrect number of material readers");
    auto ids = matids.begin();
    for(auto & mread : matreaders) {
      auto & idlist = *ids++;
      if(mread.head.first == i) {
        idlist.emplace_back(mread.head.first, mread.head.second);
        mread.advance();
      }
    }
  }

  virtual util::gid num_entities(entity_kind<D> k) const override {
    flog_assert(k == entity_kind<D>::cells || k == entity_kind<D>::vertices ||
                  k == entity_kind<D>::faces,
      "invalid entity kind " << k);

    switch(k) {
      case entity_kind<D>::cells:
        return header.elements;
        break;
      case entity_kind<D>::vertices:
        return header.nodes;
        break;
      case entity_kind<D>::faces:
        return f2v.size();
        break;
      default:
        flog_fatal("unsupported entity kind");
        break;
    }
  }

  virtual std::tuple<std::vector<util::point<D>>, std::optional<bnd_ids>>
  vertex_data(const iota_view & rng) const override {
    std::vector<util::point<D>> ret;
    ret.reserve(rng.size());

    std::optional<bnd_ids> boundaries;
    for(const auto & bread : bndreaders) {
      if(bread.head <= rng.back()) {
        if(!boundaries.has_value()) {
          boundaries.emplace();
          (*boundaries).resize(bndreaders.size());
        }
      }
    }

    for(auto i : rng) {
      if(boundaries.has_value())
        check_boundaries(*boundaries, i);

      emplace_node(ret, i, std::make_index_sequence<D>());
    }

    return std::tuple(std::move(ret), std::move(boundaries));
  }

protected:
  void check_boundaries(bnd_ids & bndids, util::gid i) const {
    flog_assert(bndids.size() == bndreaders.size(),
      "Incorrect number of boundary readers");
    auto ids = bndids.begin();
    for(auto & bread : bndreaders) {
      auto & idlist = *ids++;
      if(bread.head == i) {
        idlist.emplace_back(bread.head);
        bread.advance();
      }
    }
  }

  template<std::size_t... I>
  void emplace_node(std::vector<util::point<D>> & a,
    std::size_t i,
    std::index_sequence<I...>) const {
    util::point<D> p{nodes[i][I]...};
    a.emplace_back(p);
  }

  void read_all_nodes() {
    std::size_t pos = fh.tellg();
    flog_assert(
      pos == seeker->node_pos, "X3D: incorrect position in file for nodes.");
    std::string tok;
    nodes.resize(header.nodes);
    for(std::size_t i = 0; i < header.nodes; ++i) {
      fh >> tok; // node id
      for(std::size_t d = 0; d < D; ++d) {
        fh >> nodes[i][d];
      }
      std::getline(fh, tok); // read rest of line
    }
    std::getline(fh, tok);
    assert_str("end_nodes", tok);
  }

  x3d_face parse_face() {
    x3d_face ret;

    fh >> ret.id;
    --ret.id;
    fh >> ret.num_nodes;
    ret.nodes.resize(ret.num_nodes);
    for(std::size_t i = 0; i < ret.num_nodes; ++i) {
      fh >> ret.nodes[i];
      --ret.nodes[i];
    }
    std::string tok;
    fh >> tok; // owner processor
    fh >> tok; // neighbor processor
    fh >> ret.neighbor;
    --ret.neighbor;
    fh >> tok >> tok >> tok >> tok >> tok; // read reserved columns
    std::getline(fh, tok); // advance line

    return ret;
  }

  void read_all_faces() {
    std::size_t pos = fh.tellg();
    flog_assert(
      pos == seeker->face_pos(), "X3D: incorrect position in file for faces.");
    std::string tok;
    std::getline(fh, tok);
    assert_str("faces", tok);

    std::size_t nfaces = 0;
    for(std::size_t i = 0; i < header.faces; ++i) {
      auto face = parse_face();
      if(face.id < face.neighbor) {
        f2v.add_row(face.nodes);
        face_m2p[face.id] = nfaces++;
      }
      else {
        face_m2p[face.id] = ~(face_m2p[face.neighbor]);
      }
    }
    std::getline(fh, tok);
    assert_str("end_faces", tok);

    seeker->cell_pos = fh.tellg();
  }

  x3d_cell parse_cell() {
    x3d_cell ret;

    fh >> ret.id;
    --ret.id;
    fh >> ret.num_faces;
    ret.faces.resize(ret.num_faces);
    for(std::size_t i = 0; i < ret.num_faces; ++i) {
      size fid;
      fh >> fid;
      --fid;
      ret.faces[i] = face_m2p[fid];
    }

    return ret;
  }

  void read_all_cells() {
    std::size_t pos = fh.tellg();
    flog_assert(
      pos == seeker->cell_pos, "X3D: incorrect position in file for cells.");
    std::string tok;
    std::getline(fh, tok);
    assert_str(tok, "cells");

    for(size i = 0; i < header.elements; ++i) {
      auto cell = parse_cell();
      c2f.add_row(cell.faces);
    }

    std::getline(fh, tok); // newline
    std::getline(fh, tok);
    assert_str("end_cells", tok);

    face_m2p.clear();
  }

  x3d_header header;
  mutable std::fstream fh;
  std::optional<x3d_seek> seeker;
  std::vector<std::array<double, D>> nodes;
  util::crs c2f, f2v;
  std::unordered_map<size, size> face_m2p;
  mutable std::vector<boundary_file_reader> bndreaders;
  mutable std::vector<region_file_reader> matreaders;
};

template<std::size_t D>
std::unique_ptr<definition_base<D>>
x3d_handler(const std::string & fname,
  std::optional<std::vector<std::string>> matfiles,
  std::optional<std::vector<std::string>> bndfiles,
  MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  if(rank == 0) {
    return std::make_unique<x3d_definition<D>>(fname, matfiles, bndfiles);
  }
  else {
    return std::make_unique<undefined_definition<D>>();
  }
}
const inline bool register_x3d_2d_ =
  io_factory<2>::instance().register_type("x3d", x3d_handler<2>);
const inline bool register_x3d_3d_ =
  io_factory<3>::instance().register_type("x3d", x3d_handler<3>);

} // namespace flsp::topo::unstructured::io

#endif // FLSP_TOPO_UNSTRUCTURED_IO_X3D_DEFINITION_H
