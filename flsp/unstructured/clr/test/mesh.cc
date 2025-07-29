#include "mesh.hh"

#include "flsp/unstructured/io/exodus_definition.hh"
#include "flsp/unstructured/io/simple_definition.hh"
#include "flsp/unstructured/io/types.hh"
#include "flsp/unstructured/io/x3d_definition.hh"

#include <flecsi/util/mpi.hh>
#include <flecsi/util/unit.hh>

#include <string>
#include <vector>

using namespace unit;

template<std::size_t D>
int
mesh_initialization_test(flecsi::scheduler & sch) {
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
    typename mesh<D>::ptr m;
    sch.allocate(m,
      typename mesh<D>::mpi_coloring(
        sch, sch.runtime(), filename, matfiles, bndfiles, user_data),
      user_data);
  }; // UNIT
} // mesh_initialization_test

int
mesh_initialization(flecsi::scheduler & sch) {
  UNIT() {
    EXPECT_EQ(mesh_initialization_test<1>(sch), 0);
    EXPECT_EQ(mesh_initialization_test<2>(sch), 0);
    EXPECT_EQ(mesh_initialization_test<3>(sch), 0);
  };
} // mesh_initialization

flecsi::util::unit::driver<mesh_initialization> driver;

// I/O support.
const inline bool register_simple_1d =
  io::io_factory<policy, 1>::instance().register_type("msh",
    io::simple_handler<policy, 1>);
const inline bool register_simple_2d =
  io::io_factory<policy, 2>::instance().register_type("msh",
    io::simple_handler<policy, 2>);
const inline bool register_simple_3d =
  io::io_factory<policy, 3>::instance().register_type("msh",
    io::simple_handler<policy, 3>);

const inline bool register_x3d_1d =
  io::io_factory<policy, 2>::instance().register_type("x3d",
    io::x3d_handler<policy, 2>);
const inline bool register_x3d_2d =
  io::io_factory<policy, 3>::instance().register_type("x3d",
    io::x3d_handler<policy, 3>);

const inline bool register_exodusii_1d =
  io::io_factory<policy, 1>::instance().register_type("exo",
    io::exodus_handler<policy, 1>);
const inline bool register_exodusii_2d =
  io::io_factory<policy, 2>::instance().register_type("exo",
    io::exodus_handler<policy, 2>);
const inline bool register_exodusii_3d =
  io::io_factory<policy, 3>::instance().register_type("exo",
    io::exodus_handler<policy, 3>);
