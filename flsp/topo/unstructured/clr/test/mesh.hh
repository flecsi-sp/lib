#ifndef FLSP_TOPO_UNSTRUCTURED_CLR_TEST_MESH_HH
#define FLSP_TOPO_UNSTRUCTURED_CLR_TEST_MESH_HH

#include "flsp/topo/unstructured/clr/coloring_utils.hh"
#include "flsp/topo/unstructured/io/definition_base.hh"
#include "flsp/topo/unstructured/io/models.hh"
#include "flsp/topo/unstructured/io/types.hh"
#include "flsp/topo/unstructured/util/common.hh"

#include <flecsi/data.hh>
#include <flecsi/topo/unstructured/interface.hh>
#include <flecsi/util/parmetis.hh>

#include <cstddef>

namespace clr = flsp::topo::unstructured::clr;
namespace io = flsp::topo::unstructured::io;
namespace util = flsp::topo::unstructured::util;

/*----------------------------------------------------------------------------*
  Helper.
 *----------------------------------------------------------------------------*/

namespace unit {

using allocation = std::tuple<std::size_t, std::size_t, bool>;

template<std::size_t D>
struct policy;

template<>
struct policy<1> : flecsi::topo::help {
  static constexpr std::size_t dimension() {
    return 1;
  }
  // clang-format off
  enum index_space {
    vertices,
    edges = vertices,
    faces = edges,
    cells
  };
  static auto num_index_spaces() { return 2; }
  // clang-format on
  using index_spaces = has<vertices, cells>;
  using connectivities =
    list<from<cells, to<vertices>>, from<vertices, to<cells>>>;
  enum entity_list { owned, exclusive, shared, ghost };
  using entity_lists =
    list<entity<vertices, has<owned, exclusive, shared, ghost>>,
      entity<cells, has<owned, exclusive, shared, ghost>>>;
  static std::vector<allocation> allocations() {
    return {{index_spaces::index<vertices>, index_spaces::index<cells>, true},
      {index_spaces::index<cells>, index_spaces::index<vertices>, false}};
  }
  struct user_data {
    std::vector<util::point<1>> coords;
    std::vector<util::crs> c2v;
  };
};

template<>
struct policy<2> : flecsi::topo::help {
  static constexpr std::size_t dimension() {
    return 2;
  }
  enum index_space { vertices, edges, faces = edges, cells };
  static auto num_index_spaces() {
    return 3;
  }
  // this defines the idx_spaces order
  using index_spaces = has<vertices, edges, cells>;
  using connectivities = list<
    // this defines the cnx_allocs order
    from<vertices, to<edges, cells>>,
    from<edges, to<vertices, cells>>,
    from<cells, to<vertices, edges>>>;
  enum entity_list { owned, exclusive, shared, ghost };
  using entity_lists =
    list<entity<vertices, has<owned, exclusive, shared, ghost>>,
      entity<edges, has<owned, exclusive, shared, ghost>>,
      entity<cells, has<owned, exclusive, shared, ghost>>>;
  static std::vector<allocation> allocations() {
    return {{index_spaces::index<vertices>, index_spaces::index<edges>, true},
      {index_spaces::index<vertices>, index_spaces::index<cells>, true},
      {index_spaces::index<edges>, index_spaces::index<vertices>, false},
      {index_spaces::index<edges>, index_spaces::index<cells>, true},
      {index_spaces::index<cells>, index_spaces::index<vertices>, false},
      {index_spaces::index<cells>, index_spaces::index<edges>, false}};
  }
  struct user_data {
    std::vector<util::point<2>> coords;
    std::vector<util::crs> c2v;
    std::vector<util::crs> c2e;
    std::vector<util::crs> e2v;
  };
};

template<>
struct policy<3> : flecsi::topo::help {
  static constexpr std::size_t dimension() {
    return 3;
  }
  enum index_space { vertices, edges, faces, cells };
  static auto num_index_spaces() {
    return 4;
  }
  using index_spaces = has<vertices, edges, faces, cells>;
  using connectivities = list<from<vertices, to<edges, faces, cells>>,
    from<edges, to<vertices, cells>>,
    from<faces, to<vertices, cells>>,
    from<cells, to<vertices, edges, faces>>>;
  enum entity_list { owned, exclusive, shared, ghost };
  using entity_lists =
    list<entity<vertices, has<owned, exclusive, shared, ghost>>,
      entity<edges, has<owned, exclusive, shared, ghost>>,
      entity<cells, has<owned, exclusive, shared, ghost>>>;
  static std::vector<allocation> allocations() {
    return {{index_spaces::index<vertices>, index_spaces::index<edges>, true},
      {index_spaces::index<vertices>, index_spaces::index<faces>, true},
      {index_spaces::index<vertices>, index_spaces::index<cells>, true},
      {index_spaces::index<edges>, index_spaces::index<vertices>, false},
      {index_spaces::index<edges>, index_spaces::index<cells>, true},
      {index_spaces::index<faces>, index_spaces::index<vertices>, false},
      {index_spaces::index<faces>, index_spaces::index<cells>, true},
      {index_spaces::index<cells>, index_spaces::index<vertices>, false},
      {index_spaces::index<cells>, index_spaces::index<edges>, false},
      {index_spaces::index<cells>, index_spaces::index<faces>, false}};
  }
  struct user_data {
    std::vector<util::point<3>> coords;
    std::vector<util::crs> c2v;
    std::vector<util::crs> c2e;
    std::vector<util::crs> c2f;
    std::vector<util::crs> f2v;
    std::vector<util::crs> e2v;
  };
};

/*----------------------------------------------------------------------------*
  Mesh Specailization.
 *----------------------------------------------------------------------------*/

template<std::size_t D>
struct mesh
  : policy<D>,
    flecsi::topo::specialization<flecsi::topo::unstructured, mesh<D>> {

  static_assert(D == 1 || D == 2 || D == 3, "invalid mesh dimension");

  /*--------------------------------------------------------------------------*
    Useful Types.
   *--------------------------------------------------------------------------*/

  using typename mesh<D>::specialization::base;
  using typename mesh<D>::specialization::core;
  using coloring =
    typename flecsi::topo::specialization<flecsi::topo::unstructured,
      mesh<D>>::coloring;
  using point = util::point<D>;
  using Color = flecsi::Color;
  using entity_kind = io::entity_kind<D>;

  /*--------------------------------------------------------------------------*
    Policy Information.
   *--------------------------------------------------------------------------*/

  using index_space = typename policy<D>::index_space;
  using index_spaces = typename policy<D>::index_spaces;
  using connectivities = typename policy<D>::connectivities;
  using entity_list = typename policy<D>::entity_list;
  using entity_lists = typename policy<D>::entity_lists;
  using user_data = typename policy<D>::user_data;

  template<auto>
  static constexpr flecsi::PrivilegeCount privilege_count = 3;

  /*--------------------------------------------------------------------------*
    Expose all entity lists names directly
   *--------------------------------------------------------------------------*/

  static constexpr auto owned = entity_list::owned;
  static constexpr auto exclusive = entity_list::exclusive;
  static constexpr auto shared = entity_list::shared;
  static constexpr auto ghost = entity_list::ghost;

  /*--------------------------------------------------------------------------*
    Accessor interface.
   *--------------------------------------------------------------------------*/

  template<typename B>
  struct interface : B {
    /*!
      Return the range of all cells.

      @code
      for(auto c: m.cells()) {
        flog(trace) << f[c] << std::endl;
      } // for
      @endcode
     */
    FLECSI_INLINE_TARGET auto cells() const {
      return B::template entities<index_space::cells>();
    }

    /*!
      Return the specified cell range. Valid ranges are:
        - owned
        - exclusive
        - shared
        - ghost

      @code
      for(auto c: m.cells<owned>()) {
        f[c] += 1.0;
      } // for
      @endcode
     */
    template<typename B::entity_list L>
    FLECSI_INLINE_TARGET auto cells() const {
      return B::template special_entities<index_space::cells, L>();
    }

    /*!
      Return the range of cells from the given entity id.

      @code
      for(auto v: m.vertices<owned>()) {
        for(auto c: m.cells(v)) {
          f[v] += g[c];
        } // for
      } // for
      @endcode
     */
    template<index_space From>
    FLECSI_INLINE_TARGET auto cells(flecsi::topo::id<From> from) const {
      return B::template entities<index_space::cells>(from);
    }

    /*!
      Return the range of all vertices.

      @code
      for(auto v: m.vertices()) {
        flog(trace) << f[v] << std::endl;
      } // for
      @endcode
     */
    FLECSI_INLINE_TARGET auto vertices() const {
      return B::template entities<index_space::vertices>();
    }
    template<typename B::entity_list L>
    FLECSI_INLINE_TARGET auto vertices() const {
      return B::template special_entities<index_space::vertices, L>();
    }
    template<index_space F>
    FLECSI_INLINE_TARGET auto vertices(flecsi::topo::id<F> from) const {
      return B::template entities<index_space::vertices>(from);
    }

    FLECSI_INLINE_TARGET auto edges() const {
      return B::template entities<index_space::edges>();
    }
    template<typename B::entity_list L>
    FLECSI_INLINE_TARGET auto edges() const {
      return B::template special_entities<index_space::edges, L>();
    }
    template<index_space F>
    FLECSI_INLINE_TARGET auto edges(flecsi::topo::id<F> from) const {
      return B::template entities<index_space::edges>(from);
    }

    FLECSI_INLINE_TARGET auto faces() const {
      return B::template entities<index_space::faces>();
    }
    template<typename B::entity_list L>
    FLECSI_INLINE_TARGET auto faces() const {
      return B::template special_entities<index_space::faces, L>();
    }
    template<index_space F>
    FLECSI_INLINE_TARGET auto faces(flecsi::topo::id<F> from) const {
      return B::template entities<index_space::faces>(from);
    }

#if 0 // FIXME: Do I need these?
    template<index_space R, index_space F, index_space G>
    FLECSI_INLINE_TARGET auto intersect(flecsi::topo::id<F> from1,
      flecsi::topo::id<G> from2) const {
      using id_t = flecsi::topo::id<R>;
      auto res1 = B::template entities<R>(from1);
      auto res2 = B::template entities<R>(from2);
      std::vector<id_t> vres1(res1.begin(), res1.end());
      std::vector<id_t> vres2(res2.begin(), res2.end());
      std::sort(vres1.begin(), vres1.end());
      std::sort(vres2.begin(), vres2.end());
      std::vector<id_t> out(std::min(res1.size(), res2.size()));
      auto out_end = std::set_intersection(
        vres1.begin(), vres1.end(), vres2.begin(), vres2.end(), out.begin());
      out.resize(std::distance(out.begin(), out_end));
      return out;
    }
    template<index_space I>
    FLECSI_INLINE_TARGET auto ispace() const {
      return B::template entities<I>();
    }
    template<index_space I, index_space F>
    FLECSI_INLINE_TARGET auto ispace(flecsi::topo::id<F> from) const {
      return B::template entities<I>(from);
    }
    template<index_space I, typename B::entity_list L>
    FLECSI_INLINE_TARGET auto ispace() const {
      return B::template special_entities<I, L>();
    }
#endif
  }; // struct interface

  /*--------------------------------------------------------------------------*
    Coloring method.
   *--------------------------------------------------------------------------*/

  static coloring color(Color num_colors,
    std::string const & filename,
    std::vector<std::string> const & matfiles,
    std::vector<std::string> const & bndfiles,
    user_data & user_data) {

    std::unique_ptr<io::definition_base<D>> md =
      io::make_definition<D>(filename, matfiles, bndfiles);

    auto [rank, size] = util::mpi::info(MPI_COMM_WORLD);
    auto global_cells = util::mpi::one_to_allv(
      [&md](int, int) { return md->num_entities(entity_kind::cells); });
    auto global_vertices = util::mpi::one_to_allv(
      [&md](int, int) { return md->num_entities(entity_kind::vertices); });

    util::equal_map cem(global_cells, size);
    util::equal_map vem(global_vertices, size);
    util::equal_map pem(num_colors, size);

    /*------------------------------------------------------------------------*
      Read cell data on root process and send to the naive owners.
     *------------------------------------------------------------------------*/

    auto [c2v, finfo, minfo] = util::mpi::one_to_alli(
      [&md, &cem ](int r, int) -> auto { return md->cell_data(cem[r]); });

    /*------------------------------------------------------------------------*
      Create cell-to-cell graph.
     *------------------------------------------------------------------------*/

    auto v2c = clr::cells_through_vertices(
      cem, vem, c2v); // global vertex to global cell id
    auto [c2c, naive] = clr::create_naive(cem, c2v, v2c, D);

    /*------------------------------------------------------------------------*
      Partition with ParMETIS.
     *------------------------------------------------------------------------*/

    std::vector<Color> cell_raw;
    cell_raw = flecsi::util::parmetis::color(cem, naive, num_colors);

    /*------------------------------------------------------------------------*
      Migrate the cell data to the owning processes.
     *------------------------------------------------------------------------*/

    auto [cells, cp2m, cm2p] =
      clr::migrate_cells(cem, pem, cell_raw, c2v, finfo, minfo, c2c, v2c);

    std::map<Color, std::uint32_t> cog2l;
    std::vector<Color> col2g;

    for(auto co : pem[rank]) {
      cog2l[co] = col2g.size();
      col2g.push_back(co);
    }

    /*------------------------------------------------------------------------*
      Close the cells with respect to the halo depth.
     *------------------------------------------------------------------------*/

    coloring coloring;
    coloring.colors = num_colors;
    coloring.idx_spaces.resize(policy<D>::num_index_spaces());

    // clang-format off
    auto [vdeps, cshr, cghst, crghost, c2co, color_peers, cell_pcdata] =
      clr::close_cells(cem, pem, cell_raw, cells, 1, c2v, finfo, minfo, cp2m,
        cm2p, c2c, v2c, coloring, entity_kind::cells);
    // clang-format on

    /*------------------------------------------------------------------------*
      Assign vertex colors.
     *------------------------------------------------------------------------*/

    std::map<util::gid, Color> v2co;
    for(std::uint32_t lco{0}; lco < cells.size(); ++lco) {
      auto const & pc = coloring.idx_spaces[entity_kind::cells].colors[lco];

      for(auto c_lid /* local cell id */ : pc.owned()) {
        util::gid c = cell_pcdata[lco].all[c_lid]; // global cell id
        for(auto v : c2v[cm2p.at(c)]) {
          Color co = std::numeric_limits<Color>::max();
          for(auto cv : v2c.at(v)) {
            co = std::min(c2co.at(cv), co);
          } // for

          v2co[v] = co;
        } // for
      } // for
    } // for

    std::vector<std::vector<std::pair<util::gid, Color>>> vertices(size);
    for(auto const & [v, co] : v2co) {
      vertices[vem.bin(v)].push_back({v, co});
    } // for

    auto rank_colors =
      util::mpi::all_to_allv([&vertices](
                               int r, int) -> auto & { return vertices[r]; });

    std::vector<Color> vertex_raw;
    const auto vr = vem[rank];
    vertex_raw.resize(vr.size());
    for(auto r : rank_colors) {
      for(auto [id, co] : r) {
        vertex_raw[id - vr.front()] = co;
      } // for
    } // for

    auto [coords, binfo] = util::mpi::one_to_alli(
      [&md, &vem ](int r, int) -> auto { return md->vertex_data(vem[r]); });

    /*------------------------------------------------------------------------*
      Migrate the vertex data to the owning processes.
     *------------------------------------------------------------------------*/

    auto [vp2m, vm2p] =
      clr::migrate_vertices(vem, pem, vertex_raw, coords, binfo);

    std::vector</* over index spaces */
      std::vector</* over local colors */
        std::vector<util::crs>>> /* using local indices */
      connectivity;

    connectivity.resize(policy<D>::num_index_spaces());
    for(uint32_t from = 0; from < policy<D>::num_index_spaces(); ++from) {
      connectivity[from].resize(pem[rank].size());
      for(auto & cnx : connectivity[from]) {
        cnx.resize(policy<D>::num_index_spaces());
      } // for
    } // for

    /*------------------------------------------------------------------------*
      Close the vertices with respect to the given halo depth (mediated
      through the cells).
     *------------------------------------------------------------------------*/

    // clang-format off
    auto [vertex_pcdata] = clr::close_vertices(vem, pem, v2co, vdeps, cells,
      c2v, cm2p, vertex_raw, cell_pcdata, coords, binfo, vm2p, vp2m, coloring,
      connectivity, color_peers, entity_kind::cells, entity_kind::vertices);
    // clang-format on

    if constexpr(D != 1) {
      /*----------------------------------------------------------------------*
        Create cell region for auxiliary creation.
       *----------------------------------------------------------------------*/

      std::vector<util::gid> cfa;
      auto const & cclr = coloring.idx_spaces[entity_kind::cells];
      {
        util::id lco{0};
        for(auto const & index_coloring : cclr.colors) {
          for(auto c /* local cell id */ : index_coloring.owned()) {
            cfa.emplace_back(cell_pcdata[lco].all[c]);
          }
          ++lco;
        } // for
      }
      cfa.insert(cfa.end(), crghost.begin(), crghost.end());

      /*
        cfap2m and cfam2p are the forward and reverse maps for lc2i, i.e., they
        map between local and global offsets for cell-to-auxiliary connectivity
        information.
       */
      std::map<util::gid, util::id> cfam2p;
      std::vector<util::gid> cfap2m(cfa.size());
      util::id ci{0};
      for(auto c : cfa) {
        cfap2m[ci] = c;
        cfam2p[c] = ci++;
      } // for

      util::crs lc2i, li2d;
      std::optional<std::pair<util::crs, util::crs>> lcn;
      if(finfo.has_value()) {
        std::vector<util::gid> lrow;
        for(auto c : cfa) {
          auto const & row = finfo->c2f[cm2p.at(c)];

          std::map<util::gid, util::id> im2p;
          util::id lfid{0};
          for(auto f : finfo->p2m) {
            im2p.try_emplace(f, lfid++);
          } // for

          lrow.resize(row.size());
          lfid = 0;
          for(auto f : row) {
            util::gid fl = im2p.at(util::get_id(f));
            lrow[lfid++] = util::sign_bit(f) ? ~fl : fl;
          } // for

          lc2i.add_row(lrow);
        } // for

        std::swap(li2d, finfo->f2v);
        lcn.emplace(std::move(lc2i), std::move(li2d));
      } // if

      /*------------------------------------------------------------------------*
        Mapping variables for whatever auxiliaries are created.
       *------------------------------------------------------------------------*/

      std::map<entity_kind, util::crs> auxmap;
      std::map<entity_kind, std::vector<clr::process_color_data> &>
        other_pcdata;

      /*------------------------------------------------------------------------*
        Add interface auxiliaries.
       *------------------------------------------------------------------------*/

      // clang-format off
      auto [c2i, i2d, ia2a, il2g, ig2l, i_pcdata] =
        clr::add_auxiliaries<D, io::interface_kind<D>(), clr::heuristic::cells>(
          pem, cfa, c2v, cm2p, cfam2p, cfap2m, cshr, cghst, coloring, cog2l,
          col2g, c2co, v2co, auxmap, lcn, cell_pcdata, vertex_pcdata);

      auxmap.try_emplace(io::interface_kind<D>(), c2i);
      other_pcdata.emplace(io::interface_kind<D>(), i_pcdata);
      clr::convert_connectivity<D, io::interface_kind<D>()>(
        cell_pcdata, vertex_pcdata, other_pcdata, pem[rank].size(), c2i, i2d,
        ia2a, cm2p, cfam2p, ig2l, connectivity);
      // clang-format on

      /*------------------------------------------------------------------------*
        Add edges if they exist.
       *------------------------------------------------------------------------*/

      if constexpr(D == 3) {
        auto [c2e, e2d, ea2a, el2g, eg2l, e_pcdata] =
          clr::add_auxiliaries<D, entity_kind::edges>(pem,
            cfa,
            c2v,
            cm2p,
            cfam2p,
            cfap2m,
            cshr,
            cghst,
            coloring,
            cog2l,
            col2g,
            c2co,
            v2co,
            auxmap,
            std::nullopt,
            cell_pcdata,
            vertex_pcdata);

        auxmap.try_emplace(entity_kind::edges, c2e);
        other_pcdata.emplace(entity_kind::edges, e_pcdata);
        clr::convert_connectivity<D, entity_kind::edges>(cell_pcdata,
          vertex_pcdata,
          other_pcdata,
          pem[rank].size(),
          c2e,
          e2d,
          ea2a,
          cm2p,
          cfam2p,
          eg2l,
          connectivity);
      } // if
    } // if

    /*------------------------------------------------------------------------*
      Set external data.
     *------------------------------------------------------------------------*/

    // FIXME: Need to ask Philipp about the details of this
    // Why isn't `get_connectivity` a local method?
    user_data.coords = {coords};

    user_data.c2v = clr::get_connectivity(
      connectivity, entity_kind::cells, entity_kind::vertices);

    if constexpr(D == 2 || D == 3) {
      user_data.c2e = clr::get_connectivity(
        connectivity, entity_kind::cells, entity_kind::edges);
      user_data.e2v = clr::get_connectivity(
        connectivity, entity_kind::edges, entity_kind::vertices);
    } // if

    if constexpr(D == 3) {
      user_data.c2f = clr::get_connectivity(
        connectivity, entity_kind::cells, entity_kind::faces);
      user_data.f2v = clr::get_connectivity(
        connectivity, entity_kind::faces, entity_kind::vertices);
    } // if

    /*------------------------------------------------------------------------*
      Allocate space for connectivities.
     *------------------------------------------------------------------------*/

    for(auto [from, to, transpose] : policy<D>::allocations()) {
      for(auto & ic : coloring.idx_spaces[from].colors) {
        ic.cnx_allocs.resize(policy<D>::num_index_spaces());
      }
    } // for

    // Set the connectivity allocation sizes.
    for(auto [from, to, transpose] : policy<D>::allocations()) {
      // auto from_idx = idmap.at(from);
      // auto to_idx = idmap.at(to);
      auto from_idx = from;
      auto to_idx = to;
      for(std::size_t lco{0}; lco < cog2l.size(); ++lco) {
        coloring.idx_spaces[from_idx].colors[lco].cnx_allocs[to_idx] =
          connectivity[transpose ? to_idx : from_idx][lco]
                      [transpose ? from_idx : to_idx]
                        .values.size();
      } // for
    } // for

    coloring.color_peers.reserve(color_peers.size());
    for(const auto & s : color_peers) {
      coloring.color_peers.push_back(s.size());
    } // for

    flecsi::topo::concatenate(
      coloring.color_peers, pem.total(), MPI_COMM_WORLD);

    return coloring;
  } // color

  /*--------------------------------------------------------------------------*
    Initialization task.
   *--------------------------------------------------------------------------*/

  static void initialize(flecsi::data::topology_slot<mesh<D>> & s,
    coloring const & c,
    const policy<D>::user_data & user_data) {

    /*------------------------------------------------------------------------*
      Resize connectivity storage.
     *------------------------------------------------------------------------*/
    auto & c2v =
      s->template get_connectivity<index_space::cells, index_space::vertices>();
    auto & v2c =
      s->template get_connectivity<index_space::vertices, index_space::cells>();

    c2v(s).get_elements().resize();
    v2c(s).get_elements().resize();

    if constexpr(D == 2 || D == 3) {
      auto & c2f =
        s->template get_connectivity<index_space::cells, index_space::faces>();
      auto & f2v = s->template get_connectivity<index_space::faces,
        index_space::vertices>();

      c2f(s).get_elements().resize();
      f2v(s).get_elements().resize();
    } // if

    if constexpr(D == 3) {
      auto & c2e =
        s->template get_connectivity<index_space::cells, index_space::edges>();
      auto & e2v = s->template get_connectivity<index_space::edges,
        index_space::vertices>();

      c2e(s).get_elements().resize();
      e2v(s).get_elements().resize();
    } // if
  } // initialize
}; // struct mesh

} // namespace unit

#endif // FLSP_TOPO_UNSTRUCTURED_CLR_TEST_MESH_HH
