#include "flsp/topo/unstructured/io/definition_base.hh"
#include "flsp/topo/unstructured/io/simple_definition.hh"
//#include "flsp/topo/unstructured/io/x3d_definition.hh"
#include "flsp/topo/unstructured/io/types.hh"
#include "flsp/topo/unstructured/mesh/coloring_utils.hh"

#include <flecsi/util/mpi.hh>
#include <flecsi/util/parmetis.hh>
#include <flecsi/util/unit.hh>

#include <memory>
#include <string>

using namespace flecsi;
namespace base = flsp::topo::unstructured;
namespace io = base::io;

template<std::size_t D>
int
coloring() {
  UNIT("TASK") {
    std::unique_ptr<io::definition_base<D>> md;
    std::string mesh;

    if constexpr(D == 1) {
    }
    else if constexpr(D == 2) {
      mesh = "mesh2d-8x8.msh";
      md = io::make_definition<2>(mesh, {}, {});
    }
    else /* D == 3 */ {
      mesh = "mesh3d-8x8x8.x3d";
      std::vector<std::string> mats{"mesh3d-8x8x8.mat1", "mesh3d-8x8x8.mat2"};
      md = io::make_definition<3>(mesh, mats, {});
    }

    /*------------------------------------------------------------------------*
      Header Information.
     *------------------------------------------------------------------------*/

    auto [rank, size] = util::mpi::info(MPI_COMM_WORLD);
    auto global_cells = util::mpi::one_to_allv([&md](int, int) {
      return md->num_entities(base::config<D>::index_space::cells);
    });
    auto global_vertices = util::mpi::one_to_allv([&md](int, int) {
      return md->num_entities(base::config<D>::index_space::vertices);
    });
    const flecsi::Color colors{4};
    util::equal_map cem(global_cells, size);
    util::equal_map vem(global_vertices, size);
    util::equal_map pem(colors, size);

    std::vector</* over index spaces */
      std::vector</* over local colors */
        std::vector<util::crs>>> /* using local indices */
      connectivity;

    /*------------------------------------------------------------------------*
      Read cell data on root process and send to the naive owners.
     *------------------------------------------------------------------------*/
    auto [c2v, finfo, minfo] = util::mpi::one_to_alli(
      [&md, &cem ](int r, int) -> auto { return md->cell_data(cem[r]); });

    /*------------------------------------------------------------------------*
      Create cell-to-cell graph.
     *------------------------------------------------------------------------*/
    auto v2c = base::mesh::clr::cells_through_vertices(cem, vem, c2v);
    auto [c2c, naive] = base::mesh::clr::create_naive(cem, c2v, v2c, D);

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
    /*------------------------------------------------------------------------*
      Header Information.
     *------------------------------------------------------------------------*/
    std::stringstream ss;
    ss << "##################################\n";
    ss << "Mesh File: " << mesh << std::endl;
    ss << "Initial Partition" << std::endl;
    ss << "global cells: " << global_cells
       << " global vertices: " << global_vertices << std::endl;
    for(auto [co, cs] : cells) {
      ss << "co: " << co << "\n" << flog::container{cs} << std::endl;
    } // for
    ss << "##################################";
    flog(warn) << ss.str() << std::endl;
  };
} // coloring

util::unit::driver<coloring<2>> driver2;
