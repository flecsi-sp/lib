#ifndef FLSP_UNSTRUCTURED_IO_SIMPLE_DEFINITION_HH
#define FLSP_UNSTRUCTURED_IO_SIMPLE_DEFINITION_HH

#include "flsp/unstructured/io/definition_base.hh"

#include <cstddef>
#include <fstream>
#include <string>
#include <type_traits>
#include <vector>

namespace flsp::unstructured::io {

template<template<std::size_t> typename P,
  std::size_t D,
  typename = typename std::enable_if_t<(D == 1) || (D == 2) || (D == 3)>>
struct simple_definition : definition_base<P, D> {
  static_assert(required_keys<P, D, P<D>::cells, P<D>::vertices>::value,
    "required entity kind undefined");

  simple_definition(std::string const & filename) {
    file_.open(filename, std::ifstream::in);

    if(file_.good()) {
      std::string line;
      std::getline(file_, line);
      std::istringstream iss(line);

      // Read the number of vertices and cells
      iss >> num_vertices_ >> num_cells_;

      // Get the offset to the beginning of the vertices
      vertex_start_ = file_.tellg();

      for(size_t i(0); i < num_vertices_; ++i) {
        std::getline(file_, line);
      } // for

      cell_start_ = file_.tellg();
    }
    else {
      flog_fatal("failed opening " << filename);
    } // if
  } // simple_definition

  simple_definition(const simple_definition &) = delete;
  simple_definition & operator=(const simple_definition &) = delete;

  ~simple_definition() {
    file_.close();
  }

  util::gid num_entities(P<D>::index_space is) const override {
    flog_assert(
      is == P<D>::cells || is == P<D>::vertices, "invalid index space");
    return is == P<D>::vertices ? num_vertices_ : num_cells_;
  } // num_entities

  std::tuple<util::crs, std::optional<face_info>, std::optional<mat_ids>>
  cell_data(iota_view const & r) const override {
    skip_lines(cell_start_, *r.begin());

    std::string line;
    util::crs p2v;
    for(auto p : r) {
      (void)p;
      std::getline(file_, line);
      std::istringstream iss(line);
      p2v.add_row(std::vector<size_t>(
        std::istream_iterator<size_t>(iss), std::istream_iterator<size_t>()));
    } // for

    return std::make_tuple(p2v, std::nullopt, std::nullopt);
  } // cell_data

  std::tuple<std::vector<util::point<D>>, std::optional<bnd_ids>> vertex_data(
    iota_view const & r) const override {
    skip_lines(vertex_start_, *r.begin());

    std::string line;
    std::vector<util::point<D>> c;
    for(auto v : r) {
      (void)v;
      std::getline(file_, line);
      std::istringstream iss(line);
      util::point<D> p;
      for(std::size_t d{0}; d < D; ++d) {
        iss >> p[d];
      }
      c.push_back(p);
    } // for

    return std::tuple(std::move(c), std::nullopt);
  } // vertex_data

private:
  void skip_lines(std::iostream::pos_type start, std::size_t lines) const {
    file_.seekg(start);
    for(std::size_t l{0}; l < lines; ++l) {
      file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    } // for
  } // skip_lines

  mutable std::ifstream file_;
  std::size_t num_vertices_;
  std::size_t num_cells_;
  mutable std::iostream::pos_type vertex_start_;
  mutable std::iostream::pos_type cell_start_;
}; // simple_definition

template<template<std::size_t> typename P, std::size_t D>
std::unique_ptr<definition_base<P, D>>
simple_handler(const std::string & fname,
  std::optional<std::vector<std::string>>,
  std::optional<std::vector<std::string>>,
  MPI_Comm comm) {
  const auto [rank, size] = flecsi::util::mpi::info(comm);
  if(rank == 0) {
    return std::make_unique<simple_definition<P, D>>(fname);
  }
  else {
    return std::make_unique<undefined_definition<P, D>>();
  } // if
} // simple_handler

} // namespace flsp::unstructured::io

#endif // FLSP_UNSTRUCTURED_IO_SIMPLE_DEFINITION_HH
