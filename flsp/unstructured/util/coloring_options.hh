/*
  Copyright (c) 2021, Triad National Security, LLC
  All rights reserved
 */

#ifndef BURTON_UTIL_COLORING_OPTIONS_HH
#define BURTON_UTIL_COLORING_OPTIONS_HH

#include <ostream>

namespace flsp {
namespace util {

enum class framework_t { PARMETIS, ZOLTAN };
enum class method_t { BLOCK, RANDOM, RCB, RIB, GRAPH };
enum class graph_method_t { PARMETIS, PHG, Scotch };

inline std::ostream &
operator<<(std::ostream & out, const framework_t & value) {
  return out << [value]() {
#define PROCESS_VAL(p)                                                         \
  case(framework_t::p):                                                        \
    return #p;
    switch(value) {
      PROCESS_VAL(PARMETIS);
      PROCESS_VAL(ZOLTAN);
    }
#undef PROCESS_VAL
    return "INVALID";
  }();
}

inline std::ostream &
operator<<(std::ostream & out, const method_t & value) {
  return out << [value]() {
#define PROCESS_VAL(p)                                                         \
  case(method_t::p):                                                           \
    return #p;
    switch(value) {
      PROCESS_VAL(BLOCK);
      PROCESS_VAL(RANDOM);
      PROCESS_VAL(RCB);
      PROCESS_VAL(RIB);
      PROCESS_VAL(GRAPH);
    }
#undef PROCESS_VAL
    return "INVALID";
  }();
}

inline std::ostream &
operator<<(std::ostream & out, const graph_method_t & value) {
  return out << [value]() {
#define PROCESS_VAL(p)                                                         \
  case(graph_method_t::p):                                                     \
    return #p;
    switch(value) {
      PROCESS_VAL(PARMETIS);
      PROCESS_VAL(PHG);
      PROCESS_VAL(Scotch);
    }
#undef PROCESS_VAL
    return "INVALID";
  }();
}

struct coloring_options {
  framework_t framework = framework_t::PARMETIS;
  method_t method = method_t::GRAPH;
  graph_method_t graph_method = graph_method_t::PARMETIS;

  int debug_level = 0;
  int debug_memory = 0;
};
} // namespace util
} // namespace flsp

#endif
