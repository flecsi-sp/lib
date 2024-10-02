#include "mesh.hh"

#include "flsp/topo/unstructured/io/definition_base.hh"
#include "flsp/topo/unstructured/io/simple_definition.hh"
#include "flsp/topo/unstructured/io/x3d_definition.hh"
#include "flsp/topo/unstructured/io/types.hh"
#include "flsp/topo/unstructured/clr/coloring_utils.hh"
#include "flsp/topo/unstructured/util/common.hh"

#include <flecsi/util/mpi.hh>
#include <flecsi/util/parmetis.hh>
#include <flecsi/util/unit.hh>

#include <string>
#include <vector>

using namespace unit;

template<std::size_t D>
int
mesh_initialization_test() {
  std::string filename;
  std::vector<std::string> matfiles;
  std::vector<std::string> bndfiles;

  UNIT("TASK") {
    if constexpr(D == 1) {
      filename = "mesh1d-8.msh";
    }
    else if constexpr(D == 2) {
      filename = "mesh2d-8x8.x3d";
      matfiles.push_back("mesh2d-8x8.mat1");
      matfiles.push_back("mesh2d-8x8.mat2");
      bndfiles.push_back("mesh2d-8x8.bnd1");
      bndfiles.push_back("mesh2d-8x8.bnd2");
    }
    else /* D == 3 */ {
      filename = "mesh3d-8x8x8.exo";
    } // if

    typename mesh<D>::user_data user_data;
    typename mesh<D>::slot m;
    m.allocate(
      typename mesh<D>::mpi_coloring(
        4, filename, matfiles, bndfiles, user_data),
      user_data);
  }; // UNIT
} // mesh_initialization_test

int
mesh_initialization() {
  UNIT() {
    EXPECT_EQ(mesh_initialization_test<1>(), 0);
    EXPECT_EQ(mesh_initialization_test<2>(), 0);
    EXPECT_EQ(mesh_initialization_test<3>(), 0);
  };
} // mesh_initialization

flecsi::util::unit::driver<mesh_initialization> driver;
