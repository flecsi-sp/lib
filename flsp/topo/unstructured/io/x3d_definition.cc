#include "x3d_definition.hh"

namespace flsp::topo::unstructured::io {

void
x3d_header::read(std::fstream & fh) {
  fh.seekg(0);
  std::string tok;

  std::getline(fh, tok);
  assert_str(tok, X3D_TOKEN);

  std::getline(fh, tok);
  assert_str(tok, "header");

  fh >> tok >> process;
  fh >> tok >> numdim;
  fh >> tok >> materials;
  fh >> tok >> nodes;
  fh >> tok >> faces;
  fh >> tok >> elements;
  fh >> tok >> ghost_nodes;
  fh >> tok >> slaved_nodes;
  fh >> tok >> nodes_per_slave;
  fh >> tok >> nodes_per_face;
  fh >> tok >> faces_per_cell;
  fh >> tok >> node_data_fields;
  fh >> tok >> cell_data_fields;

  std::getline(fh, tok); // newline from previous line
  std::getline(fh, tok);
  assert_str(tok, "end_header");
}

x3d_seek::x3d_seek(std::fstream & fh, const x3d_header & h)
  : cell_pos{0}, header{h} {
  seek_nodes(fh);
}

void
x3d_seek::seek_nodes(std::fstream & fh) {
  fh.seekg(0);
  std::string tok;
  do {
    std::getline(fh, tok);
  } while(tok != "nodes");
  node_pos = fh.tellg();
}

} // namespace burton::io
