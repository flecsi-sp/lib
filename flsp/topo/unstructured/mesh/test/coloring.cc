#include "flsp/topo/unstructured/config.hh"
#include "flsp/topo/unstructured/io/exodus_definition.hh"
#include "flsp/topo/unstructured/io/simple_definition.hh"
#include "flsp/topo/unstructured/io/types.hh"
#include "flsp/topo/unstructured/io/x3d_definition.hh"
#include "flsp/topo/unstructured/mesh/coloring_utils.hh"
#include "flsp/topo/unstructured/util/common.hh"

#include <flecsi/flog.hh>
#include <flecsi/util/annotation.hh>
#include <flecsi/util/mpi.hh>
#include <flecsi/util/parmetis.hh>
#include <flecsi/util/unit.hh>

#include <map>
#include <optional>
#include <vector>

using namespace flecsi;
namespace spec = flsp::topo::unstructured;

#if 0
struct main_region : util::annotation::region<util::annotation::execution> {
  inline static const std::string name{"main"};
};
#endif

template<std::size_t D>
int
coloring_test() {
  using entity_kind = flsp::topo::unstructured::entity_kind<D>;
  UNIT("TASK") {
    /*------------------------------------------------------------------------*
      Mesh Definition.
     *------------------------------------------------------------------------*/
    std::unique_ptr<flsp::topo::unstructured::io::definition_base<D>> md;

    std::string mesh_file;
    if constexpr(D == 2) {
#if 0
      mesh_file = "mesh2d-8x8.x3d";
      std::vector<std::string> mat_files{"mesh2d-8x8.mat1", "mesh2d-8x8.mat2"};
      std::vector<std::string> bnd_files{"mesh2d-8x8.bnd1",
        "mesh2d-8x8.bnd2",
        "mesh2d-8x8.bnd3",
        "mesh2d-8x8.bnd4"};
#else
      mesh_file = "mesh2d-8x8.msh";
      std::vector<std::string> mat_files{};
      std::vector<std::string> bnd_files{};
#endif

      md = flsp::topo::unstructured::io::make_definition<D>(
        mesh_file, mat_files, bnd_files);
    }
    else if(D == 3) {
#if 0
#if 0
      mesh_file = "mesh3d-8x8x8.exo";
      std::vector<std::string> mat_files{};
      std::vector<std::string> bnd_files{};
#else
      mesh_file = "mesh3d-8x8x8.x3d";
      std::vector<std::string> mat_files{
        "mesh3d-8x8x8.mat1", "mesh3d-8x8x8.mat2"};
      std::vector<std::string> bnd_files{};
#endif
#else
#if 1
      mesh_file = "mesh3d-2x2x2.x3d";
      std::vector<std::string> mat_files{};
      std::vector<std::string> bnd_files{};
#else
      std::string mesh_file{"mesh3d-2x2x2.msh"};
      std::vector<std::string> mat_files{};
      std::vector<std::string> bnd_files{};
#endif
#endif

      md = flsp::topo::unstructured::io::make_definition<D>(
        mesh_file, mat_files, bnd_files);
    } // if

    /*------------------------------------------------------------------------*
      Header Information.
     *------------------------------------------------------------------------*/
    auto [rank, size] = util::mpi::info(MPI_COMM_WORLD);

    auto global_cells = util::mpi::one_to_allv([&md](int, int) {
      return md->num_entities(
        flsp::topo::unstructured::config<D>::index_space::cells);
    });
    auto global_vertices = util::mpi::one_to_allv([&md](int, int) {
      return md->num_entities(
        flsp::topo::unstructured::config<D>::index_space::vertices);
    });
    const flecsi::Color colors{4};
    util::equal_map cem(global_cells, size);
    util::equal_map vem(global_vertices, size);
    util::equal_map pem(colors, size);

    std::vector</* over index spaces */
      std::vector</* over local colors */
        std::vector<util::crs>>> /* using local indices */
      connectivity;

    // util::annotation::rguard<main_region> main_guard;

    /*------------------------------------------------------------------------*
      Read cell data on root process and send to the naive owners.
     *------------------------------------------------------------------------*/
    auto [c2v, finfo, minfo] = util::mpi::one_to_alli(
      [&md, &cem ](int r, int) -> auto { return md->cell_data(cem[r]); });

    /*------------------------------------------------------------------------*
      Create cell-to-cell graph.
     *------------------------------------------------------------------------*/
    auto v2c = flsp::topo::unstructured::mesh::clr::cells_through_vertices(
      cem, vem, c2v);
    auto [c2c, naive] =
      flsp::topo::unstructured::mesh::clr::create_naive(cem, c2v, v2c, D);

    std::vector<flecsi::Color> cell_raw;

    /*------------------------------------------------------------------------*
      Partition with ParMETIS.
     *------------------------------------------------------------------------*/
    if(naive.size() == 1 && naive[0].size() == 0) {
      // workaround until flecsi/flecsi!475 is available
      cell_raw = {0};
    }
    else {
      cell_raw = flecsi::util::parmetis::color(cem, naive, colors);
    } // if

    /*------------------------------------------------------------------------*
      Migrate the cell data to the owning processes.
     *------------------------------------------------------------------------*/
    auto [cells, cp2m, cm2p] =
      flsp::topo::unstructured::mesh::clr::migrate_cells(
        cem, pem, cell_raw, c2v, finfo, minfo, c2c, v2c);

    // FIXME: Remove after debug
#if 1
    std::stringstream ss;
    ss << "##################################\n";
    ss << "Mesh File: " << mesh_file << std::endl;
    ss << "Initial Partition" << std::endl;
    ss << "global cells: " << global_cells
       << " global vertices: " << global_vertices << std::endl;
    for(auto [co, cs] : cells) {
      ss << "co: " << co << "\n" << flog::container{cs} << std::endl;
    } // for
    ss << "##################################";
    flog(warn) << ss.str() << std::endl;
#endif

    std::map<Color, std::uint32_t> cog2l;
    std::vector<Color> col2g;

    for(auto const & [co, cs] : cells) {
      cog2l[co] = col2g.size();
      col2g.push_back(co);
    }

#if 0
    /*------------------------------------------------------------------------*
      Close the cells with respect to the halo depth.
     *------------------------------------------------------------------------*/
    flsp::topo::unstructured::mesh::clr::coloring coloring;
    coloring.colors = colors;
    coloring.idx_spaces.resize(mesh::num_index_spaces<D>());
    coloring.idx_spaces[mesh::kind_id<D, entity_kind::cells>()]
      .entities = global_cells;
    coloring.idx_spaces[mesh::kind_id<D, entity_kind::vertices>()]
      .entities = global_vertices;

    connectivity.resize(mesh::num_index_spaces<D>());
    for(auto from : {mesh::kind_id<D, entity_kind::cells>(),
          mesh::kind_id<D, entity_kind::vertices>()}) {
      connectivity[from].resize(pem[rank].size());
      for(auto & cnx : connectivity[from]) {
        cnx.resize(mesh::num_index_spaces<D>());
      }
    }

    // clang-format off
    auto [vdeps, cshr, cghst, crghost, c2co, color_peers, cell_pcdata] =
      mesh::clr::close_cells(cem, pem, cell_raw, cells, 1, c2v, finfo,
        minfo, cp2m, cm2p, c2c, v2c, coloring,
        mesh::kind_id<D, entity_kind::cells>());
    // clang-format on

    /*------------------------------------------------------------------------*
      Assign vertex colors.
     *------------------------------------------------------------------------*/
    std::map<util::gid, Color> v2co;
    for(std::uint32_t lco{0}; lco < cells.size(); ++lco) {
      auto const & pc =
        coloring.idx_spaces[mesh::kind_id<D, entity_kind::cells>()]
          .colors[lco];

      for(auto c : pc.owned()) {
        for(auto v : c2v[cm2p.at(c)]) {
          Color co = std::numeric_limits<Color>::max();
          for(auto cv : v2c.at(v)) {
            co = std::min(c2co.at(cv), co);
          } // for

          v2co[v] = co;
        } // for
      } // for
    } // for

    std::vector<std::vector<std::pair<util::gid, Color>>> vertices(size);
    for(auto const & [co, v] : v2co) {
      vertices[vem.bin(co)].push_back({co, v});
    } // for

    auto rank_colors =
      util::mpi::all_to_allv([&vertices](
                               int r, int) -> auto & { return vertices[r]; });

    std::vector<Color> vertex_raw;
    const auto vr = vem[rank];
    vertex_raw.resize(vr.size());
    for(auto r : rank_colors) {
      for(auto [id, co] : r) {
        vertex_raw[id - vr.front()] = co;
      } // for
    } // for

    auto [coords, binfo] = util::mpi::one_to_alli(
      [&md, &vem ](int r, int) -> auto { return md->vertex_data(vem[r]); });

    /*------------------------------------------------------------------------*
      Migrate the vertex data to the owning processes.
     *------------------------------------------------------------------------*/
    auto [vp2m, vm2p] = mesh::clr::migrate_vertices<D>(
      vem, pem, vertex_raw, coords, binfo);

    /*------------------------------------------------------------------------*
      Close the vertices with respect to the given halo depth (mediated
      through the cells).
     *------------------------------------------------------------------------*/
    // clang-format off
    auto [vertex_pcdata] = mesh::clr::close_vertices<D>(vem, pem, v2co, vdeps, cells, c2v, cm2p,
				      vertex_raw, cell_pcdata, coords, binfo, vm2p, vp2m, coloring, connectivity, color_peers,
      mesh::kind_id<D, entity_kind::cells>(),
      mesh::kind_id<D, entity_kind::vertices>());
    // clang-format on

    mesh::debug::write_coloring<D>(
      "coloring.pvtu", coloring, {cp2m, cm2p, c2v}, {vp2m, vm2p, coords});
    /*------------------------------------------------------------------------*
      Create cell region for auxiliary creation.
     *------------------------------------------------------------------------*/
    std::vector<util::gid> cfa(crghost.begin(), crghost.end());
    auto const & cclr =
      coloring.idx_spaces[mesh::kind_id<D, entity_kind::cells>()];
    for(auto const & lco : cclr.colors) {
      for(auto c : lco.owned()) {
        cfa.emplace_back(c);
      }
    } // for

    /*
      cfap2m and cfam2p are the forward and reverse maps for lc2i, i.e., they
      map between local and global offsets for cell-to-auxiliary connectivity
      information.
     */
    std::map<util::gid, util::id> cfam2p;
    std::vector<util::gid> cfap2m(cfa.size());
    util::id ci{0};
    for(auto c : cfa) {
      cfap2m[ci] = c;
      cfam2p[c] = ci++;
    } // for

    util::crs lc2i, li2d;
    std::optional<std::pair<util::crs, util::crs>> lcn;
    if(finfo.has_value()) {
      std::vector<util::gid> lrow;
      for(auto c : cfa) {
        auto const & row = finfo->c2f[cm2p.at(c)];

        std::map<util::gid, util::id> im2p;
        util::id lfid{0};
        for(auto f : finfo->p2m) {
          im2p.try_emplace(f, lfid++);
        } // for

        lrow.resize(row.size());
        lfid = 0;
        for(auto f : row) {
          util::gid fl = im2p.at(util::get_id(f));
          lrow[lfid++] = util::sign_bit(f) ? ~fl : fl;
        } // for

        lc2i.add_row(lrow);
      } // for

      std::swap(li2d, finfo->f2v);
      lcn.emplace(std::move(lc2i), std::move(li2d));
    } // if

    std::map<entity_kind, util::crs> auxmap;
#if 1
    flog(warn) << "INTERFACE" << std::endl;
    auto [c2i, i2d, ia2a, il2g, ig2l, i_pcdata] =
      mesh::clr::add_auxiliaries<D, mesh::interface_kind<D>()>(
        pem,
        cfa,
        c2v,
        cm2p,
        cfam2p,
        cfap2m,
        cshr,
        cghst,
        coloring,
        cog2l,
        col2g,
        c2co,
        v2co,
        auxmap,
        lcn,
        cell_pcdata,
        vertex_pcdata);
#endif

#if 1
    flog(warn) << "CORNERS" << std::endl;
    auto [c2cnr, cnr2d, cnra2a, cnrl2g, cnrg2l, cnr_pcdata] =
      mesh::clr::add_auxiliaries<D, entity_kind::corners>(pem,
        cfa,
        c2v,
        cm2p,
        cfam2p,
        cfap2m,
        cshr,
        cghst,
        coloring,
        cog2l,
        col2g,
        c2co,
        v2co,
        auxmap,
        std::nullopt,
        cell_pcdata,
        vertex_pcdata);

    auxmap.try_emplace(entity_kind::corners, c2cnr);
#endif

#if 1
    flog(warn) << "SIDES" << std::endl;
    auto [c2s, s2d, sa2a, sl2g, sg2l, s_pcdata] =
      mesh::clr::add_auxiliaries<D, entity_kind::sides>(pem,
        cfa,
        c2v,
        cm2p,
        cfam2p,
        cfap2m,
        cshr,
        cghst,
        coloring,
        cog2l,
        col2g,
        c2co,
        v2co,
        auxmap,
        std::nullopt,
        cell_pcdata,
        vertex_pcdata);
#endif

#if 1
    if(D == 3) {
      flog(warn) << "EDGES" << std::endl;
      auto [c2e, e2d, ea2a, el2g, eg2l, e_pcdta] =
        mesh::clr::add_auxiliaries<D, entity_kind::edges>(pem,
          cfa,
          c2v,
          cm2p,
          cfam2p,
          cfap2m,
          cshr,
          cghst,
          coloring,
          cog2l,
          col2g,
          c2co,
          v2co,
          auxmap,
          std::nullopt,
          cell_pcdata,
          vertex_pcdata);
    } // if
#endif

#if 1
    auto ispc{0};
    for(auto const & pc : coloring.idx_spaces) {
      auto prco{0};
      for(auto const & is : pc.colors) {
        flog(info) << "COLOR: " << prco << " INDEX SPACE: " << ispc << "\n"
                   << is << std::endl;
        ++prco;
      } // for
      ++ispc;
    } // for

    ispc = 0;
    for(auto const & pc : coloring.idx_spaces) {
      flog(info) << "INDEX SPACE: " << ispc << std::endl
                 << "PARTITIONS\n"
                 << flog::container{coloring.idx_spaces[ispc].partitions}
                 << std::endl;
      ++ispc;
    }
#endif
#endif
  }; // UNIT
} // core
#if 0

int
coloring() {
  UNIT() {
    EXPECT_EQ(coloring_test<1>(), 0);
    EXPECT_EQ(coloring_test<2>(), 0);
    EXPECT_EQ(coloring_test<3>(), 0);
  };
}

util::unit::driver<coloring> driver;
#endif
