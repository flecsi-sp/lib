#ifndef FLSP_TOPO_UNSTRUCTURED_IO_TYPES_HH
#define FLSP_TOPO_UNSTRUCTURED_IO_TYPES_HH

#include "flsp/topo/unstructured/io/definition_base.hh"
#include "flsp/topo/unstructured/util/factory.hh"

#include <optional>

#include <mpi.h>

namespace flsp::topo::unstructured::io {

// Define the I/O factory type.
template<std::size_t D>
using io_factory = flsp::topo::unstructured::util::factory<D /* dimensionality */,
  definition_base<D> /* return type */,
  std::string /* key type */,
  /* callback args ... */
  std::string const &,
  std::optional<std::vector<std::string>>,
  std::optional<std::vector<std::string>>,
  MPI_Comm>;

/*!
  Invoke the object factory to create a mesh definition.

  @param filename  The input file.
  @param comm      An optional MPI communicator.
 */

template<std::size_t D>
inline auto
make_definition(std::string const & filename, MPI_Comm comm = MPI_COMM_WORLD) {
  return io_factory<D>::instance().create(
    filename.substr(filename.find_last_of('.') + 1),
    filename,
    std::nullopt,
    std::nullopt,
    std::move(comm));
} // make_definition

/*!
  Invoke the object factory to create a mesh definition.

  @param filename  The input file.
  @param matfiles  Vector of material file names
  @param bndfiles  Vector of boundary file names
  @param comm      An optional MPI communicator.
 */

template<std::size_t D>
inline auto
make_definition(std::string const & filename,
  std::vector<std::string> matfiles,
  std::vector<std::string> bndfiles,
  MPI_Comm comm = MPI_COMM_WORLD) {
  return io_factory<D>::instance().create(
    filename.substr(filename.find_last_of('.') + 1),
    filename,
    std::move(matfiles),
    std::move(bndfiles),
    std::move(comm));
} // make_definition

} // namespace flsp::io

#endif // FLSP_TOPO_UNSTRUCTURED_IO_TYPES_HH
