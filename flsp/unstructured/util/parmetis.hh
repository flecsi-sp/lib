/*
  Copyright (c) 2021, Triad National Security, LLC
  All rights reserved
 */

#pragma once

#include <flecsi/util/color_map.hh>
#include <flecsi/util/crs.hh>
#include <flecsi/util/mpi.hh>
#include <flecsi/util/types.hh>

#include <parmetis.h>

#include "flsp/unstructured/util/coloring_options.hh"

#include <vector>

using parmetis_real_t = real_t;

namespace flsp {
namespace util {
namespace parmetis {
/// \addtogroup utils
/// \{

inline auto
with_zero(const flecsi::util::offsets & o) {
  std::vector<idx_t> ret;
  ret.reserve(o.size() + 1);
  ret.push_back(0);
  auto & v = o.ends();
  ret.insert(ret.end(), v.begin(), v.end());
  return ret;
}

template<typename T, typename U>
std::vector<T>
as(std::vector<U> const & v) {
  return {v.begin(), v.end()};
} // as

/// Generate a coloring of the given naive graph partition into \em colors
/// colors.  This function uses \c ParMETIS_V3_PartKway.  Each
/// process in the comm must participate.
/// \param dist distribution of entities over ranks
/// \param graph local connectivity graph
/// \param colors The number of partitions to create.
/// \param comm   An MPI_Comm object that defines the number of processes.
template<std::size_t D>
inline std::vector<flecsi::Color>
color(const flecsi::util::offsets & dist,
  const flecsi::util::crs & graph,
  flecsi::Color colors,
  const flsp::util::coloring_options & co = flsp::util::coloring_options(),
  MPI_Comm comm = MPI_COMM_WORLD) {

  auto [rank, size] = flecsi::util::mpi::info(comm);

  flog_assert(dist.size() == size_t(size),
    "distribution size (" << dist.size() << ") must equal comm size(" << size
                          << ")");

  idx_t wgtflag = 0;
  idx_t numflag = 0;
  idx_t ncon = 1;
  std::vector<parmetis_real_t> tpwgts(ncon * colors, 1.0 / colors);

  // We may need to expose some of the ParMETIS configuration options.
  std::vector<parmetis_real_t> ubvec(ncon, 1.05);
  idx_t options[3]{};
  idx_t edgecut;

  std::vector<idx_t> part(graph.size());

  std::vector<idx_t> vtxdist = with_zero(dist), xadj = with_zero(graph.offsets);
  std::vector<idx_t> adjncy = as<idx_t>(graph.values);
  // ParMETIS rejects certain trivial cases and nullptr+[0,0) ranges.
  if(adjncy.empty())
    adjncy.emplace_back();

  auto sub =
    flecsi::util::mpi::comm::split(comm, part.empty() ? MPI_UNDEFINED : 0);

  if(sub) {
    idx_t parmetis_colors = colors;
    // clang-format off
    int result = ParMETIS_V3_PartKway(&vtxdist[0], &xadj[0], &adjncy[0],
        nullptr, nullptr, &wgtflag, &numflag, &ncon, &parmetis_colors, &tpwgts[0],
        ubvec.data(), options, &edgecut, part.data(), &sub.c);
    // clang-format on

    flog_assert(result == METIS_OK, "ParMETIS_V3_PartKway returned " << result);
  }

  return {part.begin(), part.end()};
} // color

/// \}
} // namespace parmetis
} // namespace util
} // namespace flsp
