/*
  Copyright (c) 2021, Triad National Security, LLC
  All rights reserved
 */

#pragma once

#include "flsp/config.hh"

#include <cstddef>
#include <flecsi/util/color_map.hh>
#include <flecsi/util/crs.hh>
#include <flecsi/util/mpi.hh>
#include <flecsi/util/types.hh>

#include <string>
#ifdef FLECSI_SP_ENABLE_ZOLTAN
#include <zoltan_cpp.h>
#endif

#include "flsp/unstructured/util/coloring_options.hh"

#include <vector>

namespace flsp {
namespace util {
namespace zoltan {
/// \addtogroup utils
/// \{

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
  const coloring_options & co = flsp::util::coloring_options(),
  MPI_Comm comm = MPI_COMM_WORLD) {

#ifdef FLECSI_SP_ENABLE_ZOLTAN
  auto [rank, size] = flecsi::util::mpi::info(comm);

  flog_assert(dist.size() == size_t(size),
    "distribution size (" << dist.size() << ") must equal comm size(" << size
                          << ")");

  Zoltan zz(comm);

  int changes, num_gid_entries, num_lid_entries, num_import, num_export, *import_procs, *import_to_part, *export_procs, *export_to_part;
  ZOLTAN_ID_PTR import_global_ids, import_local_ids, export_global_ids, export_local_ids;

  switch(co.method) {
  case flsp::util::method_t::BLOCK:
    zz.Set_Param("LB_METHOD", "BLOCK");
    break;
  case flsp::util::method_t::RANDOM:
    zz.Set_Param("LB_METHOD", "RANDOM");
    break;
  case flsp::util::method_t::RCB:
    zz.Set_Param("LB_METHOD", "RCB");
    break;
  case flsp::util::method_t::RIB:
    zz.Set_Param("LB_METHOD", "RIB");
    break;
  case flsp::util::method_t::GRAPH:
    zz.Set_Param("LB_METHOD", "GRAPH");
    switch(co.graph_method) {
    case flsp::util::graph_method_t::PARMETIS:
      zz.Set_Param("GRAPH_PACKAGE", "ParMETIS");
      break;
    case flsp::util::graph_method_t::PHG:
      zz.Set_Param("GRAPH_PACKAGE", "PHG");
      break;
          case flsp::util::graph_method_t::Scotch:
      zz.Set_Param("GRAPH_PACKAGE", "Scotch");
      break;
    }
    break;
  }
  zz.Set_Param("DEBUG_LEVEL", std::to_string(co.debug_level).c_str());
  zz.Set_Param("DEBUG_MEMORY", std::to_string(co.debug_memory).c_str());
  zz.Set_Param("LB_APPROACH", "PARTITION");
  zz.Set_Param("RETURN_LISTS", "PART");
  zz.Set_Param("AUTO_MIGRATE", "FALSE");
  zz.Set_Param("NUM_GLOBAL_PARTS", std::to_string(colors).c_str());

  struct zoltan_user_data {
    const flecsi::util::crs * graph;
    const flecsi::util::offsets * dist;
    int rank;
  } ud{&graph, &dist, rank};

  zz.Set_Num_Obj_Fn(
    [](void * data, int * ierr) -> int {
      auto ud = static_cast<const zoltan_user_data *>(data);
      *ierr = ZOLTAN_OK;
      return (*ud->dist)[ud->rank].size();
    },
    &ud);
  zz.Set_Obj_List_Fn(
    [](void * data,
      int num_gid_entries,
      int num_lid_entries,
      ZOLTAN_ID_PTR global_ids,
      ZOLTAN_ID_PTR local_ids,
      int wgt_dim,
      float * obj_wgts,
      int * ierr) {
      auto ud = static_cast<const zoltan_user_data *>(data);
      if(wgt_dim != 0) {
        *ierr = ZOLTAN_FATAL;
        return;
      }
      auto cnt = (*ud->dist)(ud->rank);
      const auto local_size = (*ud->dist)(ud->rank + 1) - cnt;
      for(int i = 0; i < local_size; ++i) {
        global_ids[i] = cnt++;
        local_ids[i] = i;
      }
      *ierr = ZOLTAN_OK;
    },
    &ud);
  zz.Set_Num_Geom_Fn([](void * data, int * ierr) -> int {
    *ierr = ZOLTAN_OK;
    return D;
  });
  // zz.Set_Geom_Fn();
  zz.Set_Num_Edges_Fn(
    [](void * data,
      int num_gid_entries,
      int num_lid_entries,
      ZOLTAN_ID_PTR global_id,
      ZOLTAN_ID_PTR local_id,
      int * ierr) -> int {
      auto ud = static_cast<const zoltan_user_data *>(data);
      *ierr = ZOLTAN_OK;
      return ud->graph->offsets[*local_id].size();
    },
    &ud);
  zz.Set_Edge_List_Fn(
    [](void * data,
      int num_gid_entries,
      int num_lid_entries,
      ZOLTAN_ID_PTR global_id,
      ZOLTAN_ID_PTR local_id,
      ZOLTAN_ID_PTR nbor_global_id,
      int * nbor_procs,
      int wgt_dim,
      float * ewgts,
      int * ierr) {
      auto ud = static_cast<const zoltan_user_data *>(data);
      if(wgt_dim != 0) {
        *ierr = ZOLTAN_FATAL;
        return;
      }
      int i = 0;
      for(auto a : (*ud->graph)[*local_id]) {
        nbor_global_id[i] = a;
        nbor_procs[i] = ud->dist->bin(a);
        ++i;
      }
      *ierr = ZOLTAN_OK;
    },
    &ud);

  int result = zz.LB_Partition(changes,
    num_gid_entries,
    num_lid_entries,
    num_import,
    import_global_ids,
    import_local_ids,
    import_procs,
    import_to_part,
    num_export,
    export_global_ids,
    export_local_ids,
    export_procs,
    export_to_part);

  if(ZOLTAN_OK != result)
    flog_fatal("Zoltan_LB_Partition returned " << result);

  // Copy output so we can free the Zoltan arrays here.
  std::vector<flecsi::Color> part{
    &export_to_part[0], &export_to_part[num_export]};
  zz.LB_Free_Part(
    &export_global_ids, &export_local_ids, &export_procs, &export_to_part);

  return std::move(part);
#else
  flog_fatal("trying to use Zoltan domain decomposition but FleCSI-SP is not built with ENABLE_ZOLTAN");
#endif
} // color

/// \}
} // namespace zoltan
} // namespace util
} // namespace flsp
