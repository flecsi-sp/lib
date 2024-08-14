#ifndef FLSP_TOPO_UNSTRUCTURED_MESH_TYPES_HH
#define FLSP_TOPO_UNSTRUCTURED_MESH_TYPES_HH

#include "flsp/topo/unstructured/util/common.hh"

#include <flecsi/flog.hh>
#include <flecsi/topo/unstructured/types.hh>
#include <flecsi/util/color_map.hh>
#include <flecsi/util/mpi.hh>

#include <ostream>

namespace flsp::topo::unstructured::mesh {

struct coloring_definition {};

/*!
  Provide storage for the primary entity definition information that is
  populated using a mesh definition format. Primary entity definitions are
  currently always translated into entity-to-vertex relational information
  (even if this is not the native format of the mesh defition type), so
  @em e2v will always be populated. If the format support pohyhedral element
  definitions, this information will be captured in the @em e2f, @em f2v, and
  @em a_p2m variables. Otherwise, these will be empty.
 */

struct primary_definition {
  /// Primary-to-vertex relations.
  util::crs p2v;
  /// Process-to-mesh id map for primaries.
  std::vector<util::gid> prc2msh;
  /// Primary-to-face relations.
  util::crs p2f;
  /// Face-to-vertex relations.
  util::crs f2v;
  /// Process-to-mesh id map for face definitions.
  std::vector<util::gid> f_prc2msh;

  friend std::ostream &
  operator<<(std::ostream & stream, primary_definition const & p) {
    stream << "Primary-to-vertex:\n" << p.p2v << std::endl;
    stream << "Primary map (processor-to-mesh):\n"
           << flecsi::flog::container{p.prc2msh} << std::endl;
    stream << "Primary-to-face:\n" << p.p2f << std::endl;
    stream << "Face-to-vertex:\n" << p.f2v << std::endl;
    return stream;
  } // operator<<
}; // primary_definition

struct auxiliary_definition {
  /// Entity-to-definition
  util::crs e2d;
  /// Forward map: process-to-mesh
  std::vector<util::id> p2m;
  /// Reverse map: mesh-to-process
  std::unordered_map<util::gid, util::id> m2p;
}; // auxiliary_definition

} // namespace burton::mesh

#endif // FLSP_TOPO_UNSTRUCTURED_MESH_TYPES_HH
