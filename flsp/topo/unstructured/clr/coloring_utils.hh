#ifndef FLSP_TOPO_UNSTRUCTURED_CLR_COLORING_UTILS_HH
#define FLSP_TOPO_UNSTRUCTURED_CLR_COLORING_UTILS_HH

#include "flsp/topo/unstructured/clr/coloring_functors.hh"
#include "flsp/topo/unstructured/io/definition_base.hh"
#include "flsp/topo/unstructured/io/models.hh"
#include "flsp/topo/unstructured/util/common.hh"

#include <flecsi/flog.hh>
#include <flecsi/topo/unstructured/types.hh>

#include <mpi.h>

#include <algorithm>
#include <cstdint>
#include <optional>

/// \cond core
namespace flsp::topo::unstructured::clr {
/// \addtogroup mesh
/// \{

template<std::size_t D>
using entity_kind = io::entity_kind<D>;

using coloring = flecsi::topo::unstructured_base::coloring;

struct process_primary_color_data {
  std::vector<util::gid> all;
  std::map<util::gid, util::id> offsets;

  auto g2l() const {
    return [&](util::gid g) { return offsets.at(g); };
  }
};

struct process_color_data : process_primary_color_data {
  std::vector<util::gid> owned;
  std::vector<util::gid> ghost;
  std::set<util::gid> shared;
  std::unordered_map<std::size_t, std::set<std::size_t>> dependents;
};

inline void
compute_interval_sizes(coloring::index_space & clrng,
  Color colors,
  MPI_Comm comm) {
  auto & local_itvls = clrng.num_intervals;
  for(auto & c : clrng.colors)
    local_itvls.push_back(c.ghost_intervals().size());
  flecsi::topo::concatenate(local_itvls, colors, comm);
}

/*!
  Build connectivity through connectivity intersection.  Given X-to-Y
  and Y-to-Z connectivities, build X-to-Z connectivity.

  @param c2f X-to-Y connectivity
  @param f2e Y-to-Z connectivity
  @return X-to-Z connectivity (c2e)
*/
inline util::crs
intersect_connectivity(const util::crs & c2f, const util::crs & f2e) {
  util::crs c2e;
  c2e.offsets.reserve(c2f.offsets.size());
  // Note: this is a rough estimate.
  c2e.values.reserve(c2f.values.size() + f2e.values.size());

  for(const auto cell : c2f) {
    std::vector<util::gid> edges;

    // accumulate edges in cell
    for(const std::size_t face : cell) {
      for(const std::size_t ei : f2e[face]) {
        auto it = std::find(edges.begin(), edges.end(), ei);
        if(it == edges.end()) {
          edges.push_back(ei);
        }
      }
    }
    c2e.add_row(edges);
  }

  return c2e;
} // intersect_connectivity

/*!
 */
inline auto
cells_through_vertices(util::equal_map const & cem,
  util::equal_map const & vem,
  util::crs const & c2v,
  MPI_Comm comm = MPI_COMM_WORLD) {
  // This will be populated with vertex-to-cell connecivity
  // (including remotes.)
  std::map<util::gid, std::vector<util::gid>> v2c;

  auto [rank, size] = util::mpi::info(comm);

  /*
    Create a map of vertex-to-cell connectivity information from
    the initial cell distribution.
   */

  // Populate local vertex connectivity information
  const std::size_t offset = cem(rank);

  std::size_t i{0};
  for(auto c : c2v) {
    for(auto v : c) {
      v2c[v].emplace_back(offset + i);
    } // for
    ++i;
  } // for

  auto referencers =
    util::mpi::all_to_allv<vertex_referencers>({v2c, vem, rank}, comm);

  /*
    Update our local connectivity information. We now have all
    vertex-to-cell connectivity informaiton for the naive distribution
    of cells that we own.
   */

  i = 0;
  for(auto & r : referencers) {
    for(auto v : r) {
      for(auto c : v.second) {
        v2c[v.first].emplace_back(c);
      } // for
    } // for
    ++i;
  } // for

  // Remove duplicate referencers
  util::unique_each(v2c);

  /*
    Invert the vertex referencer information, i.e., find the cells
    that are referenced by our connected vertices.
   */

  std::vector<std::vector<util::gid>> referencer_inverse(size);

  for(auto const & v : v2c) {
    for(auto c : v.second) {
      const int r = cem.bin(c);
      if(r != rank) {
        referencer_inverse[r].emplace_back(v.first);
      } // if
    } // for
  } // for

  // Remove duplicate inverses
  util::unique_each(referencer_inverse);

  // Request vertex-to-cell connectivity for the cells that are
  // on other processes in the naive cell distribution.
  auto connectivity = util::mpi::all_to_allv<cell_connectivity>(
    {referencer_inverse, v2c, vem, rank}, comm);

  for(auto & r : connectivity) {
    for(auto & v : r) {
      for(auto c : v.second) {
        v2c[v.first].emplace_back(c);
      } // for
    } // for
  } // for

  // Remove duplicate referencers
  util::unique_each(v2c);

  return v2c;
} // cells_through_vertices

/*!
  Create a naive partitioning of the cells that is suitable to be passed to a
  partitioning tool. Cells are "connected" if they share \em shared vertices,
  which corresponds to the topological dimension of the connectivity, e.g,.
  cells that share at least three vertices are connected through
  two-dimensional faces.

  @param cem    An equal_map specifying the cell graph distribution.
  @param c2v    Cell-to-vertex connectivity.
  @param v2c    Vertex-to-cell connectivity.
  @param shared The number of shared vertices that constitute connectivity.
  @param comm   An MPI_Comm to use for distributed-memory communication.

  @return A std::pair<std::map<gid, std::vector<gid>>>, util::crs> object with
  the local cell-to-cell connectivity graph.
 */
inline auto
create_naive(util::equal_map const & cem,
  util::crs const & c2v,
  std::map<util::gid, std::vector<util::gid>> const & v2c,
  util::id shared,
  MPI_Comm comm = MPI_COMM_WORLD) {
  std::map<util::gid, std::vector<util::gid>> c2c;

  auto [rank, size] = util::mpi::info(comm);

  /*
    Fill in the entity-to-entity connectivity that have "shared" vertices
    in common.
   */

  std::size_t c = cem(rank);
  for(auto const & cdef : c2v) {
    std::map<util::gid, util::id> shr;

    for(auto v : cdef) {
      auto it = v2c.find(v);
      if(it != v2c.end()) {
        for(auto rc : v2c.at(v)) {
          if(rc != c)
            ++shr[rc];
        } // for
      } // if
    } // for

    for(auto tc : shr) {
      if(tc.second > shared - 1) {
        c2c[c].emplace_back(tc.first);
        c2c[tc.first].emplace_back(c);
      } // if
    } // for

    ++c;
  } // for

  // Remove duplicate connections
  util::unique_each(c2c);

  /*
    Populate the actual distributed crs data structure.
   */

  util::crs naive;

  for(const std::size_t c : cem[rank]) {
    auto it = c2c.find(c);
    if(it != c2c.end()) {
      naive.add_row(it->second);
    }
    else {
      // Unconnected cells still need a zero-length for the partitioner to work.
      naive.add_row(std::vector<util::gid>());
    }
  } // for

  return std::make_pair(std::move(c2c), std::move(naive));
} // create_graph

/*!
  Return owner information for the given @em request.

  @param em      equal_map of entity type (naive coloring).
  @param pem     equal_map of processes (colors over size).
  @param request The global ids of the desired entities.
  @param raw     The raw (un-migrated) coloring of the entities.
  @param comm    An optional MPI communicator.
 */
inline std::vector<Color>
request_owners(util::equal_map const & em,
  util::equal_map const & pem,
  std::vector<util::gid> const & request,
  std::vector<Color> const & raw,
  MPI_Comm comm = MPI_COMM_WORLD) {

  auto [rank, size] = util::mpi::info(comm);

  std::vector<std::vector<util::gid>> requests(size);

  for(auto e : request) {
    requests[em.bin(e)].emplace_back(e);
  } // for

  auto requested = util::mpi::all_to_allv(
    [&requests](int r, int) -> auto & { return requests[r]; }, comm);

  std::vector<std::vector<util::gid>> fulfill(size);
  {
    const std::size_t start = em(rank);
    Color r{0};
    for(auto rv : requested) {
      for(auto e : rv) {
        fulfill[r].emplace_back(pem.bin(raw[e - start]));
      } // for
      ++r;
    } // for
  } // scope

  auto fulfilled = util::mpi::all_to_allv(
    [&fulfill](int r, int) -> auto & { return fulfill[r]; }, comm);

  std::vector<std::size_t> offs(size, 0ul);
  std::vector<Color> owners;
  for(auto e : request) {
    auto p = em.bin(e);
    owners.emplace_back(fulfilled[p][offs[p]++]);
  } // for

  return owners;
} // request_owners

/*
  Migrate cell data to the owning colors from the initial naive owners using
  the coloring specified in @em raw.

  @param cem   Cells equal_map specifying the naive coloring.
  @param pem   Process equal_map specifying the color distribution over
               processes.
  @param raw   The coloring.
  @param c2v   Cell-to-vertex connectivity of the naive coloring. This variable
               will be modified in-place with the migrated connectivity
               information for the owned cells.
  @param finfo Optional face information from mesh definition. This variable
               will be modified in-place with migrated face information.
  @param minfo Optional material information from mesh definition. This
               variable will be modified in-place with migrated material
               information.
  @param c2c   Cell-to-cell connectivity of the naive coloring. This variable
               will be modified in-place with the migrated connectivity
               information for the owned cells.
  @param v2c   Vertex-to-cell connectivity of the naive coloring. This variable
               will be modified in-place with the migrated connectivity
               information for the owned cells.
  @param comm  An optional MPI communicator.
 */
inline auto
migrate_cells(util::equal_map const & cem,
  util::equal_map const & pem,
  std::vector<Color> const & raw,
  util::crs & c2v,
  std::optional<io::face_info> & finfo,
  std::optional<io::mat_ids> & minfo,
  std::map<util::gid, std::vector<util::gid>> & c2c,
  std::map<util::gid, std::vector<util::gid>> & v2c,
  MPI_Comm comm = MPI_COMM_WORLD) {

  auto [rank, size] = util::mpi::info(comm);

  auto migrated = util::mpi::all_to_allv<move_cells>(
    {cem, pem, raw, c2v, finfo, minfo, c2c, v2c, rank}, comm);

  // Update local information.
  std::map<util::gid, util::id> m2p;
  std::vector<util::gid> p2m;
  std::map<Color, std::vector<util::gid>> cells;
  std::set<util::gid> lf;
  for(auto const & [cell_pack, fi_pack, mi_pack, v2c_pack, c2c_pack] :
    migrated) {
    for(auto const & [info /* co, id */,
          def /* cell definition */,
          deps /* unused here */] : cell_pack) {
      auto [co, id] = info;
      c2v.add_row(def);
      m2p[id] = c2v.size() - 1;
      p2m.emplace_back(id);
      cells[co].emplace_back(id);
    } // for

    if(fi_pack.has_value()) {
      auto const & fip = *fi_pack;
      auto & fi = *finfo;

      // c2f
      for(auto const & c : fip.c2f) {
        fi.c2f.add_row(c);
      } // for

      // f2v & p2m
      auto bdef = fip.f2v.begin();
      for(auto const & f : fip.p2m) {
        auto const & def = *bdef++;

        // Unique faces
        auto [it, ret] = lf.insert(f);
        if(ret) {
          fi.p2m.emplace_back(f);
          fi.f2v.add_row(def);
        } // if
      } // for
    } // if

    if(mi_pack.has_value()) {
      auto const & mip = *mi_pack;
      auto & mi = *minfo;
      mi.resize(mip.size());
      std::size_t i{0};
      for(auto const & m : mip) {
        for(auto const & p : m) {
          mi[i].emplace_back(p);
        } // for
        ++i;
      } // for
    } // if

    for(auto const & [id, cs] : v2c_pack) {
      v2c.try_emplace(id, cs);
    } // for

    for(auto const & [id, cs] : c2c_pack) {
      c2c.try_emplace(id, cs);
    } // for
  } // for

  return std::make_tuple(std::move(cells), std::move(p2m), std::move(m2p));
} // migrate_cells

/*!
  Create the dependency closure of the local (to process) cells with respect to
  the given @em depth.

  @param cem   Cell equal-map.
  @param pem   Process equal-map.
  @param raw   Naive color assignments for cells (used to request owner
               information).
  @param cells Cell color information, i.e., a map from a global color to the
               list of local cells.
  @param depth The halo depth. This determines the depth of padding of the cell
               closure.
  @param c2v   Cell-to-vertex graph.
  @param finfo Interface information.
  @param minfo Material information.
  @param p2m   Process-to-mesh map (cells).
  @param m2p   Mesh-to-process map (cells).
  @param c2c   Cell-to-cell map.
  @param v2c   Vertex-to-cell map.
  @param comm  MPI Communicator.
 */
inline auto
close_cells(util::equal_map const & cem,
  util::equal_map const & pem,
  std::vector<Color> const & raw,
  std::map<Color, std::vector<util::gid>> const & cells,
  util::id depth,
  util::crs & c2v,
  std::optional<io::face_info> & finfo,
  std::optional<io::mat_ids> & minfo,
  std::vector<util::gid> & p2m,
  std::map<util::gid, util::id> & m2p,
  std::map<util::gid, std::vector<util::gid>> & c2c,
  std::map<util::gid, std::vector<util::gid>> & v2c,
  coloring & clrng,
  util::id cidx,
  MPI_Comm comm = MPI_COMM_WORLD) {

  std::map<util::gid, Color> c2co; /* cell-to-color */
  std::map<util::gid, std::set<Color>> dependents, dependencies, vdeps;
  std::unordered_map<Color, std::set<util::gid>> shared, ghost;
  std::vector<std::size_t> rghost;

  auto [rank, size] = util::mpi::info(comm);

  std::map<Color, std::vector<util::gid>> wkset;
  for(auto const & [co /* global color */, cs /* cells */] : cells) {
    wkset[co].reserve(cs.size());
    for(auto c : cs) {
      c2co.try_emplace(c, co);
      wkset[co].emplace_back(c);
    } // for
  } // for

  for(util::id d{0}; d < depth + 1; ++d) {
    std::vector<util::gid> layer;

    for(auto const & [co /* global color */, cs] : cells) {
      for(auto const & c : wkset.at(co)) {
        for(auto const & v : c2v[m2p.at(c)]) {
          for(auto const & cn : v2c.at(v)) {
            if(c2c.find(cn) == c2c.end()) {
              // If we don't have this cell, we need to request it.
              layer.emplace_back(cn);

              // This cell depends on the current color.
              dependencies[cn].insert(co);

              // If we're within the requested depth, this is also
              // a ghost cell for the current color.
              if(d < depth) {
                ghost[co].insert(cn);
                rghost.emplace_back(cn);
              } // if
            }
            // This cell is on the local process but not owned
            // by the current color.
            else if(d < depth && c2co.at(cn) != co) {
              // This is a shared cell.
              shared[c2co.at(cn)].insert(cn);

              // Add the current color as a dependent of this cell.
              vdeps[cn].insert(co);
              dependents[cn].insert(co);

              // This cell is a ghost for the current color.
              ghost[co].insert(cn);
            } // if
          } // for
        } // for
      } // for

      wkset.at(co).clear();
    } // for

    util::force_unique(layer);

    std::vector<std::vector<std::pair<util::gid, std::set<Color>>>> request(
      size);
    {
      auto owners = request_owners(cem, pem, layer, raw, comm);

      std::size_t cid{0};
      for(auto c : layer) {
        request[owners[cid++]].emplace_back(
          std::make_pair(c, dependencies.at(c)));
      } // for
    } // scope

    auto requested = util::mpi::all_to_allv(
      [&request](int r, int) -> auto & { return request[r]; }, comm);

    std::vector<std::vector<util::gid>> fulfill(size);
    Color r{0};
    for(auto rv : requested) {
      for(auto [c /* global cell id */, deps /* dependencies */] : rv) {
        fulfill[r].emplace_back(c);

        if(d < depth) {
          dependents[c].insert(deps.begin(), deps.end());
          shared[c2co.at(c)].insert(c);
        } // if
      } // for
      ++r;
    } // for

    // clang-format off
    auto fulfilled = util::mpi::all_to_allv<communicate_cells>(
      {fulfill, dependents, c2co, c2v, d < depth ? finfo : std::nullopt,
        d < depth ? minfo : std::nullopt, c2c, v2c, m2p}, comm);
    // clang-format on

    std::set<util::gid> lf;
    if(finfo.has_value()) {
      lf.insert(finfo->p2m.begin(), finfo->p2m.end());
    } // if

    // Update local information.
    for(auto const & [cell_pack, fi_pack, mi_pack, v2c_pack, c2c_pack] :
      fulfilled) {
      for(auto const & [info /* co, id */,
            def /* cell definition */,
            deps /* dependencies */] : cell_pack) {
        auto [co, id] = info;
        c2v.add_row(def);
        m2p[id] = c2v.size() - 1;
        p2m.emplace_back(id);
        c2co.try_emplace(id, co);

        if(d < depth && deps.has_value()) {
          vdeps[id].insert(deps->begin(), deps->end());
        } // if

        for(auto co : dependencies.at(id)) {
          wkset.at(co).emplace_back(id);
        } // for
      } // for

      if(fi_pack.has_value()) {
        auto const & fip = *fi_pack;
        auto & fi = *finfo;

        // c2f
        for(auto const & c : fip.c2f) {
          fi.c2f.add_row(c);
        } // for

        // f2v & p2m
        auto bdef = fip.f2v.begin();
        for(auto const & f : fip.p2m) {
          auto const & def = *bdef++;

          // Unique faces
          auto [it, ret] = lf.insert(f);
          if(ret) {
            fi.p2m.emplace_back(f);
            fi.f2v.add_row(def);
          } // if
        } // for
      } // if

      if(mi_pack.has_value()) {
        auto const & mip = *mi_pack;
        auto & mi = *minfo;
        mi.resize(mip.size());
        std::size_t i{0};
        for(auto const & m : mip) {
          for(auto const & p : m) {
            mi[i].emplace_back(p);
          } // for
          ++i;
        } // for
      } // if

      for(auto const & [id, cs] : v2c_pack) {
        v2c.try_emplace(id, cs);
      } // for

      for(auto const & [id, cs] : c2c_pack) {
        c2c.try_emplace(id, cs);
      } // for
    } // for
  } // for

  util::force_unique(rghost);

  // Populate coloring.
  auto & cell_color = clrng.idx_spaces[cidx];
  auto & is_peers = cell_color.peers;

  cell_color.colors.resize(pem[rank].size());

  std::vector<std::vector<util::gid>> sources(size);
  std::vector<std::vector<std::pair<Color, util::gid>>> process_ghosts(size);

  auto & partitions = clrng.idx_spaces[cidx].partitions;
  std::vector<std::set<Color>> color_peers(cells.size());

  std::vector<process_primary_color_data> cell_pcdata(pem[rank].size());

  {
    util::id lco{0};
    for(auto const & [co /* global color */, cs] : cells) {
      auto & pc = cell_pcdata[lco];
      auto & cp = color_peers[lco];
      auto & offsets = cell_pcdata[lco].offsets;
      pc.all = cs;

      for(auto c : ghost[co]) {
        const auto gco = c2co.at(c);
        const auto pr = pem.bin(gco);
        pc.all.emplace_back(c);

        sources[pr].push_back(c);
        process_ghosts[pr].emplace_back(lco, c);
        cp.insert(gco);
      }

      auto & ic = cell_color.colors[lco];
      partitions.emplace_back(pc.all.size());
      ic.entities = pc.all.size();
      {
        util::id i = 0;
        for(auto g : pc.all)
          // assign local cell IDs
          offsets[g] = i++;
        std::set<Color> peers;
        for(auto e : shared[co]) {
          const util::id lid = offsets.at(e);
          for(auto d : dependents.at(e))
            ic.peers[d].shared.insert(lid);
          peers.insert(dependents.at(e).begin(), dependents.at(e).end());
        }
        is_peers.emplace_back(peers.begin(), peers.end());
        cp.merge(peers);
      }

      ++lco;
    } // for

    /*
      Communicate local offsets for shared mesh ids.
    */

    std::vector</* by process */ std::vector<util::id>> fulfills;
    for(const auto &rv : util::mpi::all_to_allv(
          [&sources](int r, int) -> auto & { return sources[r]; }, comm)) {
      auto & f = fulfills.emplace_back();
      for(const auto id : rv)
        f.emplace_back(
          cell_pcdata[c2co.at(id) - pem[rank].front()].offsets.at(id));
    } // for

    /*
      Send/Receive the local offset information with other processes.
     */

    {
      auto pgi = process_ghosts.begin();
      for(const auto &ans : util::mpi::all_to_allv(
            [f = std::move(fulfills)](int r, int) { return std::move(f[r]); },
            comm)) {
        auto ai = ans.begin();
        for(auto & [lco, e /* global id */] : *pgi++)
          cell_color.colors[lco].peers[c2co.at(e)].ghost[*ai++] =
            cell_pcdata[lco].offsets.at(e);
      }
    }

    cell_color.entities = cem.total();

  } // scope

  flecsi::topo::concatenate(is_peers, pem.total(), comm);

  flecsi::topo::concatenate(partitions, pem.total(), comm);

  compute_interval_sizes(clrng.idx_spaces[cidx], pem.total(), comm);

  return std::make_tuple(std::move(vdeps),
    std::move(shared),
    std::move(ghost),
    std::move(rghost),
    std::move(c2co),
    std::move(color_peers),
    std::move(cell_pcdata));
} // close_cells

template<std::size_t D>
auto
migrate_vertices(util::equal_map const & vem,
  util::equal_map const & pem,
  std::vector<Color> const & raw,
  std::vector<util::point<D>> & coords,
  std::optional<io::bnd_ids> & binfo,
  MPI_Comm comm = MPI_COMM_WORLD) {

  auto [rank, size] = util::mpi::info(comm);

  auto migrated = util::mpi::all_to_allv<move_vertices<D>>(
    {vem, pem, raw, coords, binfo, rank}, comm);

  std::map<util::gid, util::id> m2p;
  std::vector<util::gid> p2m;
  for(auto const & [vertex_pack, bi_pack] : migrated) {
    for(auto const & [info, crds] : vertex_pack) {
      auto [co, id] = info;
      coords.emplace_back(crds);
      m2p[id] = coords.size() - 1;
      p2m.emplace_back(id);
    } // for

    if(bi_pack.has_value()) {
      auto const & bip = *bi_pack;
      auto bi = *binfo;
      std::size_t i{0};
      for(auto const & b : bip) {
        for(auto const & p : b) {
          bi[i].emplace_back(p);
        } // for
        ++i;
      } // for
    } // if
  } // for

  return std::make_tuple(std::move(p2m), std::move(m2p));
} // migrate_vertices

/*!
 */
template<std::size_t D>
auto
close_vertices(util::equal_map const & vem,
  util::equal_map const & pem,
  std::map<util::gid, Color> & v2co,
  std::map<util::gid, std::set<Color>> vdeps,
  std::map<Color, std::vector<util::gid>> const & cells,
  util::crs const & c2v,
  std::map<util::gid, util::id> const & cm2p,
  std::vector<Color> const & raw,
  std::vector<process_primary_color_data> const & cell_pcdata,
  std::vector<util::point<D>> & coords,
  std::optional<io::bnd_ids> & binfo,
  std::map<util::gid, util::id> & vm2p,
  std::vector<util::gid> & vp2m,
  coloring & clrng,
  std::vector<std::vector<std::vector<util::crs>>> & connectivity,
  std::vector<std::set<Color>> & color_peers,
  util::id cidx,
  util::id vidx,
  MPI_Comm comm = MPI_COMM_WORLD) {

  auto [rank, size] = util::mpi::info(comm);

  auto & vert_color = clrng.idx_spaces[vidx];
  auto & vert_partitions = vert_color.partitions;
  auto & is_peers = vert_color.peers;

  vert_color.colors.resize(pem[rank].size());

  std::vector<util::gid> remote;

  std::vector<std::vector<util::gid>> sources(size);
  std::vector<std::vector<std::pair<Color, util::gid>>> process_ghosts(size);

  std::vector<process_color_data> vertex_pcdata(pem[rank].size());

  {
    util::id lco{0};
    for(auto const & [co, cs] : cells) {
      auto primary = clrng.idx_spaces[cidx].colors[lco];
      auto & vertex_pcd = vertex_pcdata[lco];
      auto & primary_pcd = cell_pcdata[lco];
      auto & offsets = vertex_pcdata[lco].offsets;

      // Go through the shared primaries and look for ghosts. Some of these
      // may be on the local processor, i.e., we don't need to request
      // remote information about them.
      for(const auto & [d, pe] : primary.peers) {
        for(auto lid : pe.shared) {
          auto sgid = primary_pcd.all[lid];
          for(auto v : c2v[cm2p.at(sgid)]) {
            auto vit = v2co.find(v);
            flog_assert(vit != v2co.end(), "invalid vertex id");
            auto gco = vit->second;

            if(gco == co) {
              vertex_pcd.shared.insert(v);
              vertex_pcd.dependents[v].insert(d);
              vertex_pcd.owned.emplace_back(v);
            }
            else {
              vertex_pcd.ghost.emplace_back(v);
              if(gco < pem[rank].front() && gco >= *pem[rank].end()) {
                // The ghost is remote: add to remote requests.
                remote.emplace_back(v);
              } // if
            } // if
          } // for
        } // for
      } // for

      // Go through the ghost primaries and look for ghosts. Some of these
      // may also be on the local processor.
      for(auto & [c, pe] : primary.peers) {
        for(auto [rid, lid] : pe.ghost) {
          auto egid = primary_pcd.all[lid];

          for(auto v : c2v[cm2p.at(egid)]) {
            auto vit = v2co.find(v);

            // Add dependents through ghosts.
            if(vertex_pcd.shared.count(v) && vdeps.count(egid)) {
              auto deps = vdeps.at(egid);
              deps.erase(co);
              vertex_pcd.dependents[v].insert(deps.begin(), deps.end());
            } // if

            if(vit != v2co.end()) {
              auto gco = vit->second;
              if(gco != co) {
                vertex_pcd.ghost.emplace_back(v);
              } // if
            }
            else {
              // The ghost is remote: add to remote requests.
              vertex_pcd.ghost.emplace_back(v);
              remote.emplace_back(v);
            } // if
          } // for
        } // for
      }

      for(auto lid : primary.exclusive()) {
        auto egid = primary_pcd.all[lid];

        for(auto v : c2v[cm2p.at(egid)]) {
          vertex_pcd.owned.emplace_back(v);
        } // for
      } // for

      util::force_unique(vertex_pcd.owned);
      util::force_unique(vertex_pcd.ghost);

      ++lco;
    } // for
  } // scope

  util::force_unique(remote);

  {
    /*
      Request colors for remote vertices
     */

    std::vector<std::vector<util::gid>> requests(size);
    for(auto e : remote) {
      requests[vem.bin(e)].emplace_back(e);
    } // for

    auto requested = util::mpi::all_to_allv(
      [&requests](int r, int) -> auto & { return requests[r]; }, comm);

    /*
      Fulfill requests from other ranks
     */

    std::vector<std::vector<Color>> fulfill(size);
    const std::size_t start = vem(rank);
    Color r = 0;
    for(const auto & rv : requested) {
      for(auto e : rv) {
        fulfill[r].emplace_back(raw[e - start]);
      } // for
      ++r;
    } // for

    auto fulfilled = util::mpi::all_to_allv(
      [&fulfill](int r, int) -> auto & { return fulfill[r]; }, comm);

    /*
      Update our local information.
     */
    std::vector<std::size_t> offs(size);
    for(auto e : remote) {
      auto p = vem.bin(e);
      v2co.try_emplace(e, fulfilled[p][offs[p]]);
      ++offs[p];
    } // for
  }

  {
    // Build local vertex offsets
    util::id lco{0};
    for(auto const & [co, cs] : cells) {
      auto & ic = vert_color.colors[lco];
      auto & vertex_pcd = vertex_pcdata[lco];
      auto & offsets = vertex_pcd.offsets;

      // Sort ghosts by remote color to improve Legion performance.

      std::map<Color, std::vector<util::gid>> ghost_by_color;

      for(auto v : vertex_pcd.ghost) {
        ghost_by_color[v2co.at(v)].emplace_back(v);
      }

      vertex_pcd.ghost.clear();
      for(auto const & [c, v] : ghost_by_color) {
        vertex_pcd.ghost.insert(vertex_pcd.ghost.end(), v.cbegin(), v.cend());
      } // for

      for(auto v : vertex_pcd.owned) {
        vertex_pcd.all.emplace_back(v);
      }

      for(auto v : vertex_pcd.ghost) {
        vertex_pcd.all.emplace_back(v);
      }

      vert_partitions.emplace_back(vertex_pcd.all.size());
      {
        util::id i = 0;
        for(auto g : vertex_pcd.all)
          offsets[g] = i++;
      }

      lco++;
    } // for
  } // scope

  {
    util::id lco{0};
    for(auto const & [co, cs] : cells) {
      auto & ic = vert_color.colors[lco];
      auto & vertex_pcd = vertex_pcdata[lco];
      auto & cp = color_peers[lco];
      auto & o = vertex_pcd.offsets;
      std::set<Color> peers;

      ic.entities = vertex_pcd.all.size();

      for(auto v : vertex_pcd.owned) {
        const util::id lid = o.at(v);

        if(vertex_pcd.shared.size() && vertex_pcd.shared.count(v)) {
          for(auto d : vertex_pcd.dependents.at(v)) {
            ic.peers[d].shared.insert(lid);
            peers.insert(d);
          }
        }
      } // for

      for(auto v : vertex_pcd.ghost) {
        const auto gco = v2co.at(v);
        const auto pr = pem.bin(gco);

        sources[pr].emplace_back(v);
        process_ghosts[pr].emplace_back(lco, v);
        cp.insert(gco);
      } // for

      cp.insert(peers.begin(), peers.end());
      is_peers.emplace_back(peers.begin(), peers.end());

      ++lco;
    } // for
  } // scope

  /*
    Communicate local offsets for shared vertex ids.
   */

  std::vector<std::vector<util::id>> fulfills;
  for(const auto &rv : util::mpi::all_to_allv(
        [&sources](int r, int) -> auto & { return sources[r]; }, comm)) {
    auto & f = fulfills.emplace_back();
    for(const auto id : rv)
      f.emplace_back(
        vertex_pcdata[v2co.at(id) - pem[rank].front()].offsets.at(id));
  } // for

  /*
    Send/Receive the local offset information with other processes.
   */

  {
    auto pgi = process_ghosts.begin();
    for(const auto &ans : util::mpi::all_to_allv(
          [f = std::move(fulfills)](int r, int) { return std::move(f[r]); },
          comm)) {
      auto ai = ans.begin();
      for(auto & [lco, e] : *pgi++)
        vert_color.colors[lco].peers[v2co.at(e)].ghost[*ai++] =
          vertex_pcdata[lco].offsets.at(e);
    }
  }

  /*
    Finish populating the vertex index coloring.
   */
  vert_color.entities = vem.total();

  /*
    Gather the tight peer information for the vertices.
   */

  flecsi::topo::concatenate(is_peers, pem.total(), comm);

  /*
    Gather partition sizes or vertices.
   */
  flecsi::topo::concatenate(vert_partitions, pem.total(), comm);

  /*
   * Compute vertex ghost interval sizes
   */

  compute_interval_sizes(clrng.idx_spaces[vidx], pem.total(), comm);

  /*
   * Build up connectivity
   */

  for(std::size_t lco = 0; lco < pem[rank].size(); ++lco) {
    auto & vertex_pcd = vertex_pcdata[lco];
    auto & primary_pcd = cell_pcdata[lco];
    auto & crs = connectivity[cidx][lco][vidx];

    for(auto egid : primary_pcd.all) {
      crs.add_row(
        flecsi::util::transform_view(c2v[cm2p.at(egid)], vertex_pcd.g2l()));
    } // for
  } // for

  return std::make_tuple(vertex_pcdata);
} // close_vertices

template<std::size_t D, entity_kind<D> K>
constexpr std::size_t
unsorted() {
  if constexpr(D == 1) {
    return 0;
  }
  else if constexpr(D == 2) {
    if constexpr(K == entity_kind<D>::edges) {
      return 0;
    }
    else if(K == entity_kind<D>::sides) {
      /* v0, v1, unsorted(c) */
      return 1;
    }
    else if(K == entity_kind<D>::corners) {
      /* v, unsorted(c) */
      return 1;
    }
  }
  else /* D == 3 */ {
    if constexpr(K == entity_kind<D>::edges) {
      return 0;
    }
    else if(K == entity_kind<D>::faces) {
      return 0;
    }
    else if(K == entity_kind<D>::sides) {
      /* v0, v1, unsorted(f, c) */
      return 2;
    }
    else if(K == entity_kind<D>::corners) {
      /* v, unsorted(c) */
      return 1;
    } // if
  } // if
} // unsorted

template<std::size_t D, entity_kind<D> K>
constexpr bool
subcell() {
  return K == entity_kind<D>::sides || K == entity_kind<D>::corners;
} // subcell

template<std::size_t D, entity_kind<D> K>
std::string
entity_kind_name() {
  if constexpr(D == 1) {
    switch(K) {
      case entity_kind<D>::cells:
        return "cells";
      case entity_kind<D>::vertices:
        return "vertices";
      default:
        flog_fatal("invalid entity kind");
    } // switch
  }
  else if constexpr(D == 2) {
    switch(K) {
      case entity_kind<D>::cells:
        return "cells";
      case entity_kind<D>::vertices:
        return "vertices";
      case entity_kind<D>::edges:
        return "edges";
      case entity_kind<D>::sides:
        return "sides";
      case entity_kind<D>::corners:
        return "corners";
      default:
        flog_fatal("invalid entity kind");
    } // switch
  }
  else /* D == 3 */ {
    switch(K) {
      case entity_kind<D>::cells:
        return "cells";
      case entity_kind<D>::vertices:
        return "vertices";
      case entity_kind<D>::edges:
        return "edges";
      case entity_kind<D>::faces:
        return "faces";
      case entity_kind<D>::sides:
        return "sides";
      case entity_kind<D>::corners:
        return "corners";
      default:
        flog_fatal("invalid entity kind");
    } // switch
  } // if
} // entity_kind_name

/*!
  Create auxiliary entities of the specified kind @em K for the given set of
  cells.

  @tparam K The entity_kind to create.

  @param cells A cells for which a unique set of auxiliaries will be created.
  @param c2v   Cell-to-vertex connectivity for (at least) the cells specified
               in @em cells.
  @param cm2p  A mesh-to-process map for the given cells.
  @param aux   Cell-to-auxiliary information that is required to create certain
               auxiliary kinds (std::optional), e.g., sides require
               cell-to-face information in 3D.
 */
template<std::size_t D, entity_kind<D> K>
auto
create_auxiliaries(std::vector<util::gid> const & cells,
  util::crs const & c2v,
  std::map<util::gid, util::id> cm2p,
  std::map<util::gid, util::id> cfam2p,
  std::map<entity_kind<D>, util::crs> const & aux,
  const std::optional<util::crs> & i2d = std::nullopt,
  const std::optional<std::map<util::gid, util::id>> & ig2l = std::nullopt) {

  util::crs a2d;
  util::crs c2a;
  std::map<entity_kind<D>, util::crs> a2a;
  std::map<std::vector<util::gid>, util::gid> def2a;
  std::vector<util::gid> sorted;
  std::vector<util::gid> these;
  static constexpr std::size_t uv = unsorted<D, K>();

  for(auto c : cells) {
    these.clear();
    auto [ce, ca2a] = io::create_cell_entities<D, K>(
      std::make_tuple(c, cm2p.at(c), cfam2p.at(c)), c2v, aux, i2d, ig2l);

    std::size_t i{0};
    for(auto const & row : ce) {
      sorted.assign(row.begin(), row.end() - uv /* random access iterator */);
      std::sort(sorted.begin(), sorted.end());

      for(std::size_t i{row.size() - uv}; i < row.size(); ++i) {
        sorted.push_back(row[i]);
      } // for

      const auto aid = def2a.size();
      if(const auto it = def2a.try_emplace(sorted, aid); it.second) {
        // This is a new auxiliary -> add it
        a2d.add_row(row);
        these.push_back(aid);

        for(auto const & [ek, a2k] : ca2a) {
          a2a[ek].add_row(a2k[i]);
        } // for
      }
      else {
        // This is an existing auxiliary -> it was created by a different cell,
        // so its orientation is reversed.
        these.push_back(
          K == io::interface_kind<D>() ? ~it.first->second : it.first->second);
      } // if

      ++i;
    } // for

    std::sort(these.begin(), these.end(), [](util::gid a, util::gid b) {
      return util::get_id(a) < util::get_id(b);
    });
    c2a.add_row(these);
  } // for

  return std::make_tuple(std::move(c2a), std::move(a2d), std::move(a2a));
} // create_auxiliaries

enum class heuristic { vertices, cells };

/*!
  Assign colors to local auxiliaries using cell or vertex coloring information.

  @tparam K Specify the entity_kind.
  @tparam H Specify the heuristic to use to assign colors. Note this parameter
            only applies to non-subcell entities. Subcell entities always use
            the enclosing cell's color.

  @param lc2a  Local cell-to-auxiliary graph.
  @param la2d  Local auxiliary-to-definition graph.
  @param cog2l Color global-to-local map.
  @param cp2m  Cell process-to-mesh map.
  @param cm2p  Cell mesh-to-process map.
  @param cshr  Shared cells.
  @param cghst Ghost cells.
  @param clrng Coloring object.
  @param cidx  Cell index space index.
  @param c2co  Cell-to-color map.
  @param v2co  Vertex-to-color map.
 */
template<std::size_t D, entity_kind<D> K, heuristic H = heuristic::vertices>
auto
color_local_auxiliaries(util::crs const & lc2a,
  util::crs const & la2d,
  std::map<Color, util::id> const & cog2l,
  std::vector<Color> col2g,
  std::vector<util::gid> const & cp2m,
  std::map<util::gid, util::id> const & cm2p,
  std::unordered_map<Color, std::set<util::gid>> const & cshr,
  std::unordered_map<Color, std::set<util::gid>> const & cghst,
  coloring const & clrng,
  std::map<util::gid, Color> const & c2co,
  std::map<util::gid, Color> const & v2co,
  std::vector<process_primary_color_data> const & cell_pcdata) {

  std::vector<std::vector<util::gid>> a2c;
  a2c.resize(la2d.size());
  std::size_t cid{0};
  for(auto row : lc2a) {
    for(auto a : row) {
      a2c[util::get_id(a)].emplace_back(cp2m.at(cid));
    } // for
    ++cid;
  } // for

  std::map<util::gid, std::pair<Color, bool>> a2co;
  std::map<util::gid, std::pair<std::set<Color>, std::vector<util::gid>>> ghost;
  std::map<Color, std::map<util::gid, std::set<Color>>> ldependents;
  std::map<Color, std::map<util::id, std::set<Color>>> ldependencies;
  util::gid cnt{0};

  util::id lco{0};
  for(auto const & pc : clrng.idx_spaces[entity_kind<D>::cells].colors) {
    auto gco = col2g[lco];
    for(util::id c_lid{0}; c_lid < pc.entities; ++c_lid) {
      util::gid c = cell_pcdata[lco].all[c_lid];
      /*
        Loop auxiliaries and assign a color using the defined heuristic.
        Subcell entity types always use the enclosing cell color.
       */
      for(auto a : lc2a[cm2p.at(c)]) {
        auto const anc /* id without complement */ = util::get_id(a);
        bool halo{false};
        Color co = std::numeric_limits<Color>::max();
        if constexpr(subcell<D, K>()) {
          co = c2co.at(c);
        }
        else {
          if constexpr(H == heuristic::vertices) {
            for(auto v : la2d[anc]) {
              co = std::min(v2co.at(v), co);
            } // for
          }
          else if constexpr(H == heuristic::cells) {
            for(auto ci : a2c[anc]) {
              co = std::min(c2co.at(ci), co);
            } // for
          } // if
        } // if

        // Keep track of halo cells
        for(auto ci : a2c[anc]) {
          if(cshr.count(gco) && cshr.at(gco).count(ci) ||
             cghst.count(gco) && cghst.at(gco).count(ci)) {
            halo = true;
            break;
          } // if
        } // for

        auto const [it, add] = a2co.try_emplace(anc, std::make_pair(co, halo));

        if(add) {
          // This auxiliary is owned by one of this process's colors.
          if(cog2l.count(co)) {
            ++cnt;
          }
          // This is a ghost auxiliary.
          else {
            auto def = util::to_vector(la2d[anc]);
            ghost[anc] = std::make_pair(std::set<Color>{gco}, std::move(def));
          } // if
        }
        else {
          if(halo) {
            a2co[it->first].second = true;
          }

          // This is a dependent.
          if(!cog2l.count(co)) {
            ghost.at(anc).first.insert(gco);
          } // if
        }

        // Keep track of ghost and shared between local colors (on process).
        if(cog2l.count(co) && co != gco) {
          // co -> owning color
          // anc -> auxiliary id (without complement)
          // pc.color -> current cell color (our color)
          ldependents[co][anc].insert(gco);
          ldependencies[co][anc].insert(gco);
        } // if
      } // for
    } // for
    ++lco;
  } // for

  return std::make_tuple(cnt,
    std::move(a2co),
    std::move(ghost),
    std::move(ldependents),
    std::move(ldependencies));
} // color_local_auxiliaries

inline auto
assign_global_ids(util::gid cnt, MPI_Comm comm = MPI_COMM_WORLD) {
  auto [rank, size] = util::mpi::info(comm);

  util::gid offset{0};
  util::gid entities{0};
  int r{0};
  for(auto c : util::mpi::all_gatherv(cnt)) {
    entities += c;
    if(r++ < rank) {
      offset += c;
    } // if
  } // for

  return std::make_pair(entities, offset);
} // assign_global_ids

/*!
  Assign colors and global ids to the auxiliaries with index @em aidx.
 */
template<std::size_t D, entity_kind<D> K>
auto
color_auxiliaries(util::equal_map const & pem,
  util::gid cnt,
  util::crs const & lc2a,
  util::crs const & la2d,
  std::map<entity_kind<D>, util::crs> const & la2a,
  std::map<util::gid,
    std::pair<std::set<Color>, std::vector<util::gid>>> const & ghost,
  std::map<Color, std::map<util::gid, std::set<Color>>> & ldependents,
  std::map<Color, std::map<util::id, std::set<Color>>> & ldependencies,
  std::map<Color, util::id> const & cog2l,
  std::map<util::gid, std::pair<Color, bool>> const & a2co,
  bool interface = false,
  MPI_Comm comm = MPI_COMM_WORLD) {

  auto [rank, size] = util::mpi::info(comm);

  // Get the total number of entities and the local offset.
  auto [entities, offset] = assign_global_ids(cnt, comm);

  std::vector<std::map<std::vector<util::gid>, std::pair<util::id, util::gid>>>
    shared(cog2l.size());
  std::vector<util::gid> sorted;
  static constexpr std::size_t uv = unsorted<D, K>();
  std::vector<util::gid> l2g(
    a2co.size(), std::numeric_limits<util::gid>::max());
  for(auto && [lid, info] : a2co) {
    auto && [co, halo] = info;

    if(cog2l.count(co)) {
      if(halo) {
        const auto & row = la2d[lid];
        sorted.assign(row.begin(), row.end() - uv);
        std::sort(sorted.begin(), sorted.end());

        for(util::gid i{row.size() - uv}; i < row.size(); ++i) {
          sorted.push_back(row[i]);
        } // for

        // Keep track of shared auxiliaries so that we can enforce consistent
        // definition ordering of ghost auxiliaires on another process. The
        // sorted definition is used on the remote to check against their local
        // ordering.
        shared[cog2l.at(co)].try_emplace(
          std::move(sorted), std::make_pair(lid, offset));
      } // if

      // This actually assigns the global color.
      l2g[lid] = offset++;
    } // if
  } // for

  std::map<Color, std::map<util::gid, std::set<Color>>> dependents;
  for(auto const & [oco, deps] : ldependents) {
    for(auto const & [a, dcos] : deps) {
      dependents[oco][l2g.at(a)].insert(dcos.begin(), dcos.end());
    } // for
  } // for
  ldependents.clear();

  std::map<Color, std::map<util::gid, std::set<Color>>> dependencies;
  for(auto const & [oco, deps] : ldependencies) {
    for(auto const & [a, dcos] : deps) {
      dependencies[oco][l2g.at(a)].insert(dcos.begin(), dcos.end());
    } // for
  } // for
  ldependencies.clear();

  std::vector<std::vector<util::id>> lids(size);
  std::vector<std::vector<std::tuple<Color /* owning color */,
    std::set<Color> /* requesting colors */,
    std::vector<util::gid>>>>
    request(size);
  for(auto const & [a, info] : ghost) {
    auto const & [dcos, def] = info;
    auto const pr{pem.bin(a2co.at(a).first)};
    request[pr].emplace_back(std::make_tuple(a2co.at(a).first, dcos, def));
    lids[pr].emplace_back(a);
  } // for

  auto requested = util::mpi::all_to_allv(
    [&request](int r, int) -> auto & { return request[r]; }, comm);

  std::vector<std::vector<std::tuple<util::gid, std::vector<util::gid>>>>
    fulfill(size);
  util::id pr{0};
  for(auto & rv : requested) {
    for(auto & [oco, dcos, def] : rv) {
      sorted.assign(def.begin(), def.end() - uv);
      std::sort(sorted.begin(), sorted.end());
      for(util::gid i{def.size() - uv}; i < def.size(); ++i) {
        sorted.push_back(def[i]);
      } // for

      auto it = shared[cog2l.at(oco)].find(sorted);
      flog_assert(it != shared[cog2l.at(oco)].end(),
        "invalid auxiliary definition (entity_kind: "
          << (entity_kind_name<D, K>()) << ")\n"
          << flog::container{sorted});

      fulfill[pr].emplace_back(std::make_tuple(
        it->second.second, util::to_vector(la2d[it->second.first])));
      dependents[oco][it->second.second].insert(dcos.begin(), dcos.end());
    } // for

    ++pr;
  } // for

  auto fulfilled = util::mpi::all_to_allv(
    [&fulfill](int r, int) -> auto & { return fulfill[r]; }, comm);

  pr = 0;
  std::map<util::gid, std::vector<util::gid>> swap;
  for(auto & fv : fulfilled) {
    std::size_t off{0};

    for(auto const & [gid, def] : fv) {
      auto & [oco, rco, rdef] = request[pr][off];
      if(def != rdef) {
        swap.try_emplace(gid, def);
      } // if

      l2g[lids[pr][off++]] = gid;
    } // for
    ++pr;
  } // for

  std::map<util::gid, util::id> g2l;
  {
    util::gid a{0};
    for(auto gid : l2g) {
      g2l.try_emplace(gid, a++);
    } // for
  } // scope

  /*
    This assigns the local ordering of the auxiliaries consistent with al2g.
   */
  // FIXME: Update these with non-const row logic when flecsi #335 is done
  util::crs a2d;
  {
    util::gid ai{0};
    for(auto const & def : la2d) {
      // The owning process has a different ordering -> use it!
      if(auto it = swap.find(l2g.at(ai)); it != swap.end()) {
        a2d.add_row(it->second);
      }
      // Our ordering is consistent or we own this auxiliary.
      else {
        a2d.add_row(def);
      } // if

      ++ai;
    } // for
  } // scope

  // FIXME: Update these with non-const value logic when flecsi #335 is done
  util::crs c2a;
  for(auto const & def : lc2a) {
    std::vector<util::gid> row(def.size());
    util::id ai{0};
    for(auto a : def) {
      decltype(a) ag = l2g.at(util::get_id(a));
      if(auto it = swap.find(ag);
         interface && (it != swap.end() || util::sign_bit(a))) {
        ag = ~ag;
      }

      row[ai] = ag;
      ++ai;
    } // for

    c2a.add_row(row);
  } // for

  for(auto const & [a, info] : ghost) {
    auto const & [dcos, def] = info;
    dependencies[a2co.at(a).first][l2g[a]].insert(dcos.begin(), dcos.end());
  } // for

  return std::make_tuple(entities,
    std::move(c2a),
    std::move(a2d),
    std::move(la2a), // currently nothing done to modify this.
    std::move(l2g),
    std::move(g2l),
    std::move(dependents),
    std::move(dependencies));
} // color_auxiliaries

inline auto
close_auxiliaries(util::equal_map const & pem,
  std::map<Color, util::id> const & cog2l,
  util::gid num_auxiliaries,
  std::map<util::gid, std::pair<Color, bool>> const & a2co,
  std::vector<util::gid> const & al2g,
  std::map<util::gid, util::id> const & ag2l,
  std::map<Color, std::map<util::gid, std::set<Color>>> const & dependents,
  std::map<Color, std::map<util::gid, std::set<Color>>> const & dependencies,
  coloring & clrng,
  std::vector<process_primary_color_data> const & cell_pcdata,
  std::vector<process_color_data> const & vertex_pcdata,
  util::id aidx,
  bool reorder = true,
  MPI_Comm comm = MPI_COMM_WORLD) {
  auto [rank, size] = util::mpi::info(comm);

  auto & aux_color = clrng.idx_spaces[aidx];

  aux_color.colors.resize(cog2l.size());

  std::vector<std::vector<util::gid>> sources(size);
  std::vector<std::vector<std::pair<Color, util::gid>>> process_ghosts(size);

  std::vector<std::set<util::gid>> pclo(cog2l.size());
  std::vector<std::set<Color>> peers(cog2l.size());

  std::vector<std::set<Color>> color_peers_(cog2l.size());

  // if(auxs_.size()) {
  //   for(std::size_t lco = 0; lco < cog2l.size(); ++lco) {
  //     pclo[lco].insert(
  //       primary_pcdata[lco].all.begin(), primary_pcdata[lco].all.end());
  //   }
  // }

  std::vector<process_color_data> aux_pcdata(cog2l.size());

  for(auto const & [gco, lco] : cog2l) {
    auto & aux_pcd = aux_pcdata[lco];
    auto & primary_pcd = cell_pcdata[lco];
    auto & vertex_pcd = vertex_pcdata[lco];
    // auto & cnx = connectivity(idx)[lco];
    auto & cp = color_peers_[lco];
    auto & offsets = aux_pcd.offsets;

    for(auto [gid, lid] : ag2l) {
      auto const co = a2co.at(lid).first;

      if(gco == co) {
        aux_pcd.owned.emplace_back(gid);

        // cnx[cell_index()].add_row(
        //   util::transform_view(aux.i2e[lid], primary_pcd.g2l()));
        // cnx[vertex_index()].add_row(
        //   util::transform_view(aux.i2v[lid], vertex_pcd.g2l()));
        // for(auto & im : cd_.aidxs) {
        //   auto it = aux.i2a.find(im.kind);
        //   if(it != aux.i2a.end()) {
        //     cnx[im.idx].add_row(util::transform_view(it->second[lid],
        //     auxiliary_pcdata(im.kind)[lco].g2l()));
        //   }
        // }

        if(!dependents.empty() && dependents.at(co).count(gid)) {
          aux_pcd.shared.insert(gid);
          auto const & deps = dependents.at(co).at(gid);

          for(auto d : deps) {
            peers[lco].insert(d);
            aux_pcd.dependents[gid].insert(d);
          }
        }
      }
      else if(!dependencies.empty() && dependencies.at(co).count(gid) &&
              dependencies.at(co).at(gid).count(gco)) {
        aux_pcd.ghost.push_back(gid);
        const auto pr = pem.bin(co);

        sources[pr].emplace_back(gid);
        process_ghosts[pr].emplace_back(lco, gid);
        cp.insert(co);

        // Only add auxiliary connectivity that is covered by
        // the primary closure.
        // std::vector<util::gid> pall;
        // for(auto e : aux.i2e[lid]) {
        //   if(pclo[lco].count(e)) {
        //     pall.push_back(e);
        //   } // if
        // } // for

        // if(pall.size()) {
        //   cnx[cell_index()].add_row(
        //                             flecsi::util::transform_view(pall,
        //                             primary_pcd.g2l()));
        // }

        // cnx[vertex_index()].add_row(
        //                             flecsi::util::transform_view(aux.i2v[lid],
        //                             vertex_pcd.g2l()));

        // for(auto & im : cd_.aidxs) {
        //   auto it = aux.i2a.find(im.kind);
        //   if(it != aux.i2a.end()) {
        //     cnx[im.idx].add_row(flecsi::util::transform_view(it->second[lid],
        //     auxiliary_pcdata(im.kind)[lco].g2l()));
        //   }
        // }
      }
    }

    util::force_unique(aux_pcd.owned);
    util::force_unique(aux_pcd.ghost);

    if(reorder) {
      // Sort ghosts by remote color to improve Legion performance.

      std::map<Color, std::vector<util::gid>> ghost_by_color;

      for(auto a : aux_pcd.ghost) {
        ghost_by_color[a2co.at(ag2l.at(a)).first].emplace_back(a);
      }

      for(auto gid : aux_pcd.owned) {
        aux_pcd.all.emplace_back(gid);
      }

      for(auto const & [c, a] : ghost_by_color) {
        aux_pcd.all.insert(aux_pcd.all.end(), a.cbegin(), a.cend());
      }
    }
    else {
      for(auto gid : aux_pcd.owned) {
        aux_pcd.all.emplace_back(gid);
      }

      for(auto gid : aux_pcd.ghost) {
        aux_pcd.all.emplace_back(gid);
      }

      util::force_unique(aux_pcd.all);
    }

    {
      util::id i = 0;
      for(auto gid : aux_pcd.all) {
        offsets[gid] = i++;
      }
    }

  } // for

  /*
    Communicate local offsets for shared aux ids.
   */

  std::vector<std::vector<util::id>> fulfills;
  for(const auto &rv : util::mpi::all_to_allv(
        [&sources](int r, int) -> auto & { return sources[r]; }, comm)) {
    auto & f = fulfills.emplace_back();
    for(const auto id : rv)
      f.emplace_back(
        aux_pcdata[cog2l.at(a2co.at(ag2l.at(id)).first)].offsets.at(id));
  } // for

  /*
    Send/Receive the local offset information with other processes.
   */

  {
    auto pgi = process_ghosts.begin();
    for(const auto &ans : util::mpi::all_to_allv(
          [f = std::move(fulfills)](int r, int) { return std::move(f[r]); },
          comm)) {
      auto ai = ans.begin();
      for(auto & [lco, e] : *pgi++)
        aux_color.colors[lco].peers[a2co.at(ag2l.at(e)).first].ghost[*ai++] =
          aux_pcdata[lco].offsets.at(e);
    }
  }

  /*
    Populate the coloring information.
   */
  aux_color.entities = num_auxiliaries;

  for(auto const & [gco, lco] : cog2l) {
    auto & ic = aux_color.colors[lco];
    auto & aux_pcd = aux_pcdata[lco];
    auto & offsets = aux_pcd.offsets;

    ic.entities = aux_pcd.all.size();

    for(auto gid : aux_pcd.owned) {
      util::id lid = offsets.at(gid);

      if(aux_pcd.shared.count(gid)) {
        for(auto d : aux_pcd.dependents[gid]) {
          ic.peers[d].shared.insert(lid);
        }
      }
    }
  } // for

  /*
    Populate peer information.
   */

  auto & parts = clrng.idx_spaces[aidx].partitions;
  auto & is_peers = clrng.idx_spaces[aidx].peers;
  for(std::size_t lco{0}; lco < cog2l.size(); ++lco) {
    auto & ic = aux_color.colors[lco];
    parts.emplace_back(ic.entities);
    auto & pp = peers[lco];
    color_peers_[lco].insert(pp.begin(), pp.end());
    is_peers.emplace_back(pp.begin(), pp.end());
  } // for

  /*
    Gather the tight peer information for this auxiliary entity type.
   */

  flecsi::topo::concatenate(is_peers, pem.total(), comm);

  /*
    Gather partition sizes for entities.
   */

  flecsi::topo::concatenate(parts, pem.total(), comm);

  /*
   * Compute aux ghost interval sizes
   */

  compute_interval_sizes(clrng.idx_spaces[aidx], pem.total(), comm);

#if 0
  for(auto const & [co, lco] : cog2l) {
    auto & ac = clrng.idx_spaces[aidx].colors[lco];
    ac.entities = num_auxiliaries;
  } // for

  std::map<Color, std::set<Color>> peers;
  for(auto gid : al2g) {
    auto const [aco, ha] = a2co.at(ag2l.at(gid));

    for(auto [co, lco] : cog2l) {
      auto & ac = clrng.idx_spaces[aidx].colors[lco];

      bool const ownd = (co == aco);
      bool const ghst = !ownd && dependencies.at(aco).count(gid) &&
                        dependencies.at(aco).at(gid).count(co);

      if(ownd) {
        ac.coloring.owned.emplace_back(gid);

        if(dependents.count(aco) && dependents.at(aco).count(gid)) {
          auto const & deps = dependents.at(co).at(gid);
          ac.coloring.shared.push_back({gid, {deps.begin(), deps.end()}});
        }
        else {
          ac.coloring.exclusive.emplace_back(gid);
        } // if
      } // if

      if(ghst) {
        auto const [pr, li] = pem.invert(aco);
        ac.coloring.ghost.push_back({gid, pr, Color(li), aco});

        peers[lco].insert(aco);
      } // if
    } // for
  } // for

  std::vector<std::size_t> & partitions = clrng.partitions[aidx];
  std::vector<std::vector<Color>> is_peers(cog2l.size());
  for(auto [co, lco] : cog2l) {
    auto & ac = clrng.idx_spaces[aidx][lco];
    partitions.emplace_back(ac.coloring.all.size());

    if(peers.count(lco)) {
      is_peers[lco].resize(peers.at(lco).size());
      std::copy(
        peers.at(lco).begin(), peers.at(lco).end(), is_peers[lco].begin());
      ac.peers.resize(peers.at(lco).size());
      std::copy(peers.at(lco).begin(), peers.at(lco).end(), ac.peers.begin());
    }
  } // for

  flecsi::topo::concatenate(partitions, pem.total(), comm);
#endif
  return std::move(aux_pcdata);
} // close_auxiliaries

/*!
  Add auxiliary with entity_kind \em K.

  @tparam K Specify the entity_kind.
  @tparam H Specify the heuristic to use to assign colors. Note this parameter
            only applies to non-subcell entities. Subcell entities always use
            the enclosing cell's color.

  @param pem    Process equal-map.
  @param cfa    Cell-for-auxilary (Footprint of cells for which to create
                auxiliaries.)
  @param c2v    Cell-to-vertex graph.
  @param cm2p   Cell mesh-to-process map.
  @param cfam2p Cell-for-auxiliary mesh-to-process map.
  @param cfap2m Cell-for-auxiliary process-to-mesh map.
  @param cshr   Shared cells.
  @param cghst  Ghost cells.
  @param clrng  Coloring object.
  @param cog2l  Color global-to-local map.
  @param c2co   Cell-to-color map.
  @param v2co   Vertex-to-color map.
  @param aux    Auxiliary kind that this auxiliary depends upon.
  @param lcn    Auxiliary information provided by the mesh format.
 */
template<std::size_t D, entity_kind<D> K, heuristic H = heuristic::vertices>
inline auto
add_auxiliaries(util::equal_map const & pem,
  std::vector<util::gid> const & cfa,
  util::crs const & c2v,
  std::map<util::gid, util::id> const & cm2p,
  std::map<util::gid, util::id> const & cfam2p,
  std::vector<util::gid> const & cfap2m,
  std::unordered_map<Color, std::set<util::gid>> const & cshr,
  std::unordered_map<Color, std::set<util::gid>> const & cghst,
  coloring & clrng,
  std::map<Color, util::id> const & cog2l,
  std::vector<Color> const & col2g,
  std::map<util::gid, Color> const & c2co,
  std::map<util::gid, Color> const & v2co,
  std::map<entity_kind<D>, util::crs> const & aux,
  std::optional<std::pair<util::crs, util::crs>> const & lcn,
  std::vector<process_primary_color_data> const & cell_pcdata,
  std::vector<process_color_data> const & vertex_pcdata,
  const std::optional<util::crs> & i2d = std::nullopt,
  const std::optional<std::map<util::gid, util::id>> & ig2l = std::nullopt) {

  // clang-format off
  util::crs lc2a, la2d;
  std::map<entity_kind<D>, util::crs> la2a;
  if(lcn.has_value()) {
    std::tie(lc2a, la2d) = std::move(*lcn);
  }
  else {
    std::tie(lc2a, la2d, la2a) =
      create_auxiliaries<D, K>(cfa, c2v, cm2p, cfam2p, aux, i2d, ig2l);
  } // if

  auto [acnt, a2co, aghost, aldependents, aldependencies] =
    color_local_auxiliaries<D, K, H>(lc2a, la2d, cog2l, col2g, cfap2m, cfam2p,
      cshr, cghst, clrng, c2co, v2co, cell_pcdata);

  auto [num_aux, c2a, a2d, a2a, al2g, ag2l, adependents, adependencies] =
    color_auxiliaries<D, K>(pem, acnt, lc2a, la2d, la2a, aghost, aldependents,
      aldependencies, cog2l, a2co, true);

  constexpr bool reorder =
    K != entity_kind<D>::faces &&
      K != entity_kind<D>::corners && K != entity_kind<D>::sides;
  
  auto aux_pcdata = close_auxiliaries(pem, cog2l, num_aux, a2co, al2g, ag2l,
    adependents, adependencies, clrng, cell_pcdata, vertex_pcdata, K, reorder);
  // clang-format on

  return std::make_tuple(std::move(c2a),
    std::move(a2d),
    std::move(a2a),
    std::move(al2g),
    std::move(ag2l),
    std::move(aux_pcdata));
} // add_auxiliaries

inline std::vector<util::crs>
get_connectivity(
  std::vector<std::vector<std::vector<util::crs>>> const & connectivity,
  uint32_t from,
  uint32_t to) {
  // FIXME: check if idmap_ should be used like in FleCSI unit test
  /* auto from_idx = idmap_.at(from); */
  /* auto to_idx = idmap_.at(to); */
  auto from_idx = from;
  auto to_idx = to;
  std::vector<util::crs> cnxs;
  for(auto & cnx : connectivity[from_idx]) {
    cnxs.push_back(cnx[to_idx]);
  }
  return cnxs;
}

template<std::size_t D, entity_kind<D> aux_kind>
inline void
convert_connectivity(
  std::vector<process_primary_color_data> const & cell_pcdata,
  std::vector<process_color_data> const & vertex_pcdata,
  std::map<entity_kind<D>, std::vector<process_color_data> &> const &
    other_pcdata,
  Color nlco,
  util::crs const & c2a,
  util::crs const & a2d,
  std::map<entity_kind<D>, util::crs> const & a2a,
  std::map<util::gid, util::id> const & cm2p,
  std::map<util::gid, util::id> const & cfam2p,
  std::map<util::gid, util::id> const & ag2l,
  std::vector<std::vector<std::vector<util::crs>>> & connectivity) {

  auto const & aux_pcdata = other_pcdata.at(aux_kind);

  for(std::size_t lco = 0; lco < nlco; ++lco) {
    auto const & primary_pcd = cell_pcdata[lco];
    auto const & vertex_pcd = vertex_pcdata[lco];
    auto const & aux_pcd = aux_pcdata[lco];
    auto const & face_pcd = other_pcdata.at(entity_kind<D>::faces)[lco];
    auto const & corner_pcd = [&]() -> decltype(auto) {
      if constexpr(aux_kind == entity_kind<D>::sides)
        return other_pcdata.at(entity_kind<D>::corners)[lco];
      else
        return process_color_data();
    }();

    auto & crs_cell = connectivity[entity_kind<D>::cells][lco][aux_kind];
    auto & crs_vert = connectivity[aux_kind][lco][entity_kind<D>::vertices];
    auto & crs_faces = connectivity[aux_kind][lco][entity_kind<D>::faces];
    auto & crs_corners = connectivity[aux_kind][lco][entity_kind<D>::corners];

    for(auto egid : primary_pcd.all) {
      crs_cell.add_row(flecsi::util::transform_view(c2a[cfam2p.at(egid)],
        [&](util::gid g) { return aux_pcd.g2l()(util::get_id(g)); }));
    } // for

    for(auto agid : aux_pcd.all) {
      const util::id alid = ag2l.at(agid);
      std::vector<util::id> vertices;
      std::vector<util::id> faces, corners;
      auto view = flecsi::util::transform_view(
        a2d[alid], [&](util::gid g) { return vertex_pcd.g2l()(g); });

      if constexpr(aux_kind == entity_kind<D>::sides) {
        auto face_view =
          flecsi::util::transform_view(a2a.at(entity_kind<D>::faces)[alid],
            [&](util::gid g) { return face_pcd.g2l()(g); });
        // auto corner_view = flecsi::util::transform_view(
        //   a2a.at(entity_kind<D>::corners)[ag2l.at(agid)],
        //   [&](util::gid g) { return corner_pcd.g2l()(g); });

        vertices.emplace_back(view[0]);
        vertices.emplace_back(view[1]);
        faces.emplace_back(face_view[0]);
        // corners.emplace_back(corner_view[0]);
        // corners.emplace_back(corner_view[1]);

        // build side-corner connectivity
        // const auto & c2cor =
        // connectivity[entity_kind::cells][lco][entity_kind::corners]; const
        // auto & cor2v =
        // connectivity[entity_kind::corners][lco][entity_kind::vertices];
        // util::id cid{0};
        // for(const auto sides : crs_cell) {
        //   if(std::find(sides.begin(), sides.end(), alid) != sides.end())
        //   break;
        //   ++cid;
        // }
        // flog_assert(cid < crs_cell.size(), "could not find cell associated
        // with side"); for(const auto cor : c2cor[cid]) {
        //   const auto vert_of_cor = cor2v[cor][0];
        //   for(const auto v : vertices) {
        //     if(v == vert_of_cor) {
        //       corners.emplace_back(cor);
        //       break;
        //     }
        //   }
        // }
        // flog_assert(corners.size() == 2, "every side needs to be connected to
        // 2 corners");
      }
      else if constexpr(aux_kind == entity_kind<D>::corners) {
        auto face_view = flecsi::util::transform_view(
          a2a.at(entity_kind<D>::faces)[ag2l.at(agid)],
          [&](util::gid g) { return face_pcd.g2l()(util::get_id(g)); });
        vertices.emplace_back(view[0]);
        for(auto f : face_view) {
          faces.emplace_back(f);
        }
      }
      else {
        vertices.assign(view.begin(), view.end());
      }
      crs_vert.add_row(vertices);
      crs_faces.add_row(faces);
      // if constexpr(aux_kind != entity_kind<D>::corners)
      //   crs_corners.add_row(corners);
    }
  }
}

/// \}
} // namespace flsp::topo::unstructured::clr
/// \endcond

#endif // FLSP_TOPO_UNSTRUCTURED_CLR_COLORING_UTILS_HH
