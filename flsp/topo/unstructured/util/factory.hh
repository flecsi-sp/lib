#ifndef FLSP_TOPO_UNSTRUCTURED_UTIL_FACTORY_HH
#define FLSP_TOPO_UNSTRUCTURED_UTIL_FACTORY_HH

#include <flecsi/flog.hh>

#include <map>
#include <memory>

namespace flsp::topo::unstructured::util {

template<std::size_t D, typename R, typename K, typename... As>
class factory
{
public:
  using handler = std::unique_ptr<R> (*)(As... args);
  using map_t = std::map<K, handler>;

  factory(const factory & f) = delete;
  factory & operator=(const factory &) = delete;

  static factory & instance() {
    static factory f;
    return f;
  } // instance

  bool register_type(const K key, const handler ch) {
    return map_.try_emplace(key, ch).second;
  } // register_type

  auto create(const K key, As &&... args) {
    typename map_t::const_iterator ita = map_.find(key);
    flog_assert(ita != map_.end(), "error unknown type " << key);
    return (ita->second(std::forward<As>(args)...));
  } // create

private:
  map_t map_;

  factory() {}
  ~factory() {}
}; // class factory

} // namespace flsp::topo::unstructured::util

#endif // FLSP_TOPO_UNSTRUCTURED_UTIL_FACTORY_HH
