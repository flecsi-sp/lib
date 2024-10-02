#ifndef FLSP_TOPO_UNSTRUCTURED_CLR_CORE_FUNCTORS_HH
#define FLSP_TOPO_UNSTRUCTURED_CLR_CORE_FUNCTORS_HH

#include "flsp/topo/unstructured/io/definition_base.hh"
#include "flsp/topo/unstructured/util/common.hh"

/// \cond core
namespace flsp::topo::unstructured::clr {
/// \addtogroup mesh
/// \{

/*----------------------------------------------------------------------------*
  Pack partial vertex-to-cell connectivity information for the rank that owns
  the vertex in the naive vertex distribution.
 *----------------------------------------------------------------------------*/
struct vertex_referencers {
  vertex_referencers(std::map<util::gid, std::vector<util::gid>> const & v2c,
    const util::equal_map & vem,
    int rank)
    : size_(vem.size()) {
    refs_.resize(size_);
    for(auto v : v2c) {
      auto r = vem.bin(v.first);
      if(int(r) != rank) {
        for(auto c : v.second) {
          refs_[r][v.first].emplace_back(c);
        } // for
      } // if
    } // for
  } // vertex_refernces

  auto & operator()(int rank, int) const {
    flog_assert(rank < size_, "invalid rank");
    return refs_[rank];
  } // operator(int, int)

private:
  const int size_;
  std::vector<std::map<util::gid, std::vector<util::gid>>> refs_;
}; // struct vertex_referencers

/*-----------------------------------------------------------------------------*
  Pack vertex-to-cell connectivity information.
 *----------------------------------------------------------------------------*/
struct cell_connectivity {
  cell_connectivity(std::vector<std::vector<util::gid>> const & p2c,
    std::map<util::gid, std::vector<util::gid>> const & v2c,
    const util::equal_map & vem,
    int rank)
    : size_(vem.size()), conn_(size_) {

    int po{0};
    for(auto & p : p2c) {
      if(po != rank) {
        for(auto v : p) {
          conn_[po][v] = v2c.at(v);
        } // for
      } // for
      ++po;
    } // for
  } // cell_connectivity

  auto & operator()(int rank, int) const {
    flog_assert(rank < size_, "invalid rank");
    return conn_[rank];
  } // operator(int, int)

private:
  const int size_;
  std::vector<std::map<util::gid, std::vector<util::gid>>> conn_;
}; // struct cell_connectivity

/*-----------------------------------------------------------------------------*
  Cell data type used for moving and communicating cell data.
 *----------------------------------------------------------------------------*/
// clang-format off
using cell_data =
  std::tuple<
    std::vector</* over cells */
      std::tuple<
        std::pair<Color, util::gid>, /* owning color and mesh id */
        std::vector<util::gid>, /* cell definition (vertex mesh ids) */
        std::optional<std::set<Color>> /* dependents */
      >
    >,
    std::optional<io::face_info>,
    std::optional<io::mat_ids>,
    std::map</* vertex-to-cell connectivity map */
      util::gid,
      std::vector<util::gid>
    >,
    std::map</* cell-to-cell connectivity map */
      util::gid,
      std::vector<util::gid>
    >
  >;
// clang-format on

/*-----------------------------------------------------------------------------*
  Move cell data to a new color.
 *----------------------------------------------------------------------------*/
struct move_cells {
  move_cells(util::equal_map const & cem,
    util::equal_map const & pem,
    std::vector<Color> const & raw,
    util::crs & c2v,
    std::optional<io::face_info> & finfo,
    std::optional<io::mat_ids> & minfo,
    std::map<util::gid, std::vector<util::gid>> & c2c,
    std::map<util::gid, std::vector<util::gid>> & v2c,
    int rank)
    : size_(cem.size()) {

    // Reverse face map for populating forward maps that we pass to other ranks.
    std::map<util::gid, util::id> f_m2p;
    if(finfo.has_value()) {
      for(std::size_t i{0}; i < finfo->p2m.size(); ++i) {
        f_m2p.try_emplace(util::get_id(finfo->p2m[i]), i);
      } // if
    } // for

    for(std::uint32_t r{0}; r < std::uint32_t(size_); ++r) {
      std::set<util::gid> lf;
      std::vector<std::tuple<std::pair<Color, util::gid>,
        std::vector<util::gid>,
        std::optional<std::set<Color>>>>
        cell_pack;
      std::optional<io::face_info> fi_pack;
      std::optional<io::mat_ids> mi_pack;
      std::map<util::gid, std::vector<util::gid>> c2c_pack;
      std::map<util::gid, std::vector<util::gid>> v2c_pack;

      if(finfo.has_value()) {
        fi_pack.emplace(io::face_info{});
      } // if

      if(minfo.has_value()) {
        mi_pack.emplace(io::mat_ids{});
        mi_pack->resize(minfo->size());
      }

      for(std::size_t i{0}; i < raw.size(); ++i) {
        if(pem.bin(raw[i] /* the color of ith local cell */) == r) {
          const auto j = cem(rank) + i;
          const std::pair<Color, util::gid> info{raw[i], j};

          // Cell information.
          cell_pack.push_back(
            std::make_tuple(info, to_vector(c2v[i]), std::nullopt));

          // Collect face information to send.
          if(finfo.has_value()) {
            fi_pack->c2f.add_row(finfo->c2f[i]);
            for(auto f : finfo->c2f[i]) {
              lf.insert(util::get_id(f));
            } // for
          } // if

          // Collect material information to send.
          if(minfo.has_value()) {
            util::gid ii{0};
            for(auto m : minfo.value()) {
              // This assumes that m is sorted.
              const auto pp = std::partition_point(m.begin(),
                m.end(),
                [j](auto const & k) { return k.first < j; });
              const auto off = std::distance(m.begin(), pp);
              if(pp != m.end() && m[off].first == j) {
                (*mi_pack)[ii].push_back(m[off]);
              } // if
              ++ii;
            } // for
          } // if

          /*
            If we have full connectivity information, we pack it up
            and send it. We do not send information that is potentially,
            or actually incomplete because additional communication
            will be required to resolve it regardless.
           */
          for(auto const & v : c2v[i]) {
            v2c_pack[v] = v2c[v];
          } // for

          c2c_pack[j] = c2c[j];

          /*
            Remove information that we are migrating. We can't remove
            vertex-to-cell information until the loop over ranks is done.
           */
          c2c.erase(i); // erase element i (not key j)
        } // if
      } // for

      // Collect the aggregate face information over the faces for this r.
      for(auto f : lf) {
        fi_pack->f2v.add_row(finfo->f2v[f_m2p.at(f)]);
        fi_pack->p2m.push_back(f);
      } // for

      packs_.emplace_back(
        std::make_tuple(cell_pack, fi_pack, mi_pack, v2c_pack, c2c_pack));
    } // for

    c2v.clear();
    c2c.clear();
    v2c.clear();

    if(finfo.has_value()) {
      finfo->c2f.clear();
      finfo->f2v.clear();
      finfo->p2m.clear();
    } // if

    if(minfo.has_value()) {
      // Don't clear the outer size => num materials.
      for(auto m : minfo.value()) {
        m.clear();
      } // for
    } // if
  } // move_cells

  const cell_data & operator()(int rank, int) const {
    flog_assert(rank < size_, "invalid rank");
    return packs_[rank];
  }

private:
  const int size_;
  std::vector<cell_data> packs_;
}; // struct move_cells

/*-----------------------------------------------------------------------------*
  Communicate cell data to a new color.
 *----------------------------------------------------------------------------*/
struct communicate_cells {
  communicate_cells(std::vector<std::vector<util::gid>> const & cells,
    std::map<util::gid, std::set<Color>> const & deps,
    std::map<util::gid, Color> const & c2co,
    util::crs const & c2v,
    std::optional<io::face_info> const & finfo,
    std::optional<io::mat_ids> const & minfo,
    std::map<util::gid, std::vector<util::gid>> const & c2c,
    std::map<util::gid, std::vector<util::gid>> const & v2c,
    std::map<util::gid, util::id> const & m2p)
    : size_(cells.size()) {

    // Reverse map for populating forward maps that we pass to other ranks.
    std::map<util::gid, util::id> f_m2p;
    if(finfo.has_value()) {
      for(std::size_t i{0}; i < finfo->p2m.size(); ++i) {
        f_m2p.try_emplace(util::get_id(finfo->p2m[i]), i);
      } // if
    } // for

    for(auto const & rv : cells) {
      std::set<util::gid> lf;
      std::vector<std::tuple<std::pair<Color, util::gid>,
        std::vector<util::gid>,
        std::optional<std::set<Color>>>>
        cell_pack;
      std::optional<io::face_info> fi_pack;
      std::optional<io::mat_ids> mi_pack;
      std::map<util::gid, std::vector<util::gid>> c2c_pack;
      std::map<util::gid, std::vector<util::gid>> v2c_pack;

      if(finfo.has_value()) {
        fi_pack.emplace(io::face_info{});
      } // if

      if(minfo.has_value()) {
        mi_pack.emplace(io::mat_ids{});
        mi_pack->resize(minfo->size());
      }

      for(auto const & c : rv) {
        const std::pair<Color, util::gid> info{c2co.at(c), c};

        auto optdep = [&deps](auto c) {
          return deps.count(c) ? std::optional<std::set<Color>>{deps.at(c)}
                               : std::nullopt;
        };

        cell_pack.push_back(
          std::make_tuple(info, to_vector(c2v[m2p.at(c)]), optdep(c)));

        // Collect face information to send.
        if(finfo.has_value()) {
          fi_pack->c2f.add_row(finfo->c2f[m2p.at(c)]);
          for(auto f : finfo->c2f[m2p.at(c)]) {
            lf.insert(util::get_id(f));
          } // for
        } // if

        // Collect material information to send.
        if(minfo.has_value()) {
          util::gid ii{0};
          for(auto m : minfo.value()) {
            // This assumes that m is sorted.
            const auto pp = std::partition_point(
              m.begin(), m.end(), [c](auto const & k) { return k.first < c; });
            const auto off = std::distance(m.begin(), pp);
            if(pp != m.end() && m[off].first == c) {
              (*mi_pack)[ii].push_back(m[off]);
            } // if
            ++ii;
          } // for
        } // if

        /*
          If we have full connectivity information, we pack it up
          and send it. We do not send information that is potentially,
          or actually incomplete because additional communication
          will be required to resolve it regardless.
         */
        for(auto const & v : c2v[m2p.at(c)]) {
          v2c_pack[v] = v2c.at(v);
        } // for

        c2c_pack[c] = c2c.at(c);
      } // for

      for(auto f : lf) {
        fi_pack->f2v.add_row(finfo->f2v[f_m2p.at(f)]);
        fi_pack->p2m.push_back(f);
      } // for

      packs_.emplace_back(
        std::make_tuple(cell_pack, fi_pack, mi_pack, v2c_pack, c2c_pack));
    } // for
  } // communicate_cells

  const cell_data & operator()(int rank, int) const {
    flog_assert(rank < size_, "invalid rank");
    return packs_[rank];
  }

private:
  const int size_;
  std::vector<cell_data> packs_;
}; // struct communicate_cells

/*-----------------------------------------------------------------------------*
  Vertex data type used for moving vertices.
 *----------------------------------------------------------------------------*/
// clang-format off
template<std::size_t D>
using vertex_data =
  std::tuple<
    std::vector<
      std::tuple<
        std::pair<Color, util::gid>,
        util::point<D>
      >
    >,
    std::optional<io::bnd_ids>
  >;
// clang-format on

/*-----------------------------------------------------------------------------*
  Move vertex data to a new color.
 *----------------------------------------------------------------------------*/
template<std::size_t D>
struct move_vertices {
  move_vertices(util::equal_map const & vem,
    util::equal_map const & pem,
    std::vector<Color> const & raw,
    std::vector<util::point<D>> & coords,
    std::optional<io::bnd_ids> const & binfo,
    int rank)
    : size_(vem.size()) {
    for(std::uint32_t r{0}; r < std::uint32_t(size_); ++r) {
      std::vector<std::tuple<std::pair<Color, util::gid>, util::point<D>>>
        vertex_pack;
      std::optional<io::bnd_ids> bi_pack;

      if(binfo.has_value()) {
        bi_pack.emplace(io::bnd_ids{});
        bi_pack->resize(binfo->size());
      } // if

      for(std::size_t i{0}; i < raw.size(); ++i) {
        if(pem.bin(raw[i]) == r) {
          const auto j = vem(rank) + i;
          const std::pair<Color, util::gid> info{raw[i], j};

          // Vertex information.
          vertex_pack.push_back(std::make_tuple(info, coords[i]));

          // Collect boundary information to send.
          if(binfo.has_value()) {
            util::gid ii{0};
            for(auto b : *binfo) {
              const auto pp = std::partition_point(
                b.begin(), b.end(), [j](auto const & k) { return k < j; });
              const auto off = std::distance(b.begin(), pp);
              if(pp != b.end() && b[off] == j) {
                (*bi_pack)[ii].push_back(j);
              } // if
              ++ii;
            } // for
          } // if
        } // if
      } // for

      packs_.emplace_back(vertex_pack, bi_pack);
    } // for

    coords.clear();

    if(binfo.has_value()) {
      // Don't clear the outer size => num boundaries.
      for(auto b : binfo.value()) {
        b.clear();
      } // for
    } // if
  } // move_vertices

  const vertex_data<D> & operator()(int rank, int) const {
    flog_assert(rank < size_, "invalid rank");
    return packs_[rank];
  } // operator()

private:
  const int size_;
  std::vector<vertex_data<D>> packs_;
}; // struct move_vertices

/*-----------------------------------------------------------------------------*
  Communicate vertex data to a new color.
 *----------------------------------------------------------------------------*/
template<std::size_t D>
struct communicate_vertices {
  communicate_vertices(std::vector<std::vector<util::gid>> const & vertices,
    std::map<util::gid, Color> const & v2co,
    std::vector<util::point<D>> const & coords,
    std::optional<io::bnd_ids> const & binfo,
    std::map<util::gid, util::id> const & m2p)
    : size_(vertices.size()) {

    for(auto const & rv : vertices) {
      std::vector<std::tuple<std::pair<Color, util::gid>, util::point<D>>>
        vertex_pack;
      std::optional<io::bnd_ids> bi_pack;

      if(binfo.has_value()) {
        bi_pack.emplace(io::bnd_ids{});
        bi_pack->resize(binfo->size());
      } // if

      for(auto const & v : rv) {
        const std::pair<Color, util::gid> info{v2co.at(v), v};
        vertex_pack.push_back(std::make_tuple(info, coords[m2p.at(v)]));

        if(binfo.has_value()) {
          util::gid ii{0};
          for(auto b : *binfo) {
            const auto pp = std::partition_point(
              b.begin(), b.end(), [v](auto const & k) { return k < v; });
            const auto off = std::distance(b.begin(), pp);
            if(pp != b.end() && b[off] == v) {
              (*bi_pack)[ii].push_back(v);
            } // if
            ++ii;
          } // for
        } // if
      } // for

      auto optbi = [&binfo, &bi_pack]() {
        return binfo.has_value() ? std::optional<io::bnd_ids>{bi_pack}
                                 : std::nullopt;
      };

      packs_.emplace_back(vertex_pack, optbi());
    } // for
  } // communicate_cells

  const vertex_data<D> & operator()(int rank, int) const {
    flog_assert(rank < size_, "invalid rank");
    return packs_[rank];
  }

private:
  const int size_;
  std::vector<vertex_data<D>> packs_;
}; // struct communicate_vertices

template<typename T>
struct move_field {

  move_field(std::size_t size,
    Color colors,
    const std::vector<Color> & index_colors,
    const std::vector<T> & field)
    : em_(colors, size), index_colors_(index_colors), field_(field) {}

  std::vector<T> operator()(int rank, int) const {
    std::vector<T> field_pack;
    int i = 0;
    for(const auto & v : index_colors_) {
      if(em_.bin(v) == static_cast<unsigned int>(rank)) {
        field_pack.push_back(field_[i]);
      } // if
      ++i;
    } // for
    return field_pack;
  } // operator(int, int)

private:
  const util::equal_map em_;
  const std::vector<Color> & index_colors_;
  const std::vector<T> & field_;
}; // struct move_field

template<typename T>
struct pack_field {
  pack_field(const util::equal_map & dist, const std::vector<T> & field)
    : dist_(dist), field_(field) {}

  auto operator()(int rank, int) const {
    std::vector<std::pair<std::size_t, T>> field_pack;
    const auto dr = dist_[rank];
    field_pack.reserve(dr.size());
    for(const std::size_t i : dr) {
      field_pack.emplace_back(i, field_[i]);
    } // for
    return field_pack;
  } // operator(int, int)

private:
  const util::equal_map & dist_;
  const std::vector<T> & field_;
}; // struct pack_field

/// \}
} // namespace flsp::topo::unstructured::clr
/// \endcond

#endif // FLSP_TOPO_UNSTRUCTURED_CLR_CORE_FUNCTORS_HH
