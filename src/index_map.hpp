#ifndef SRC_INDEX_MAP_HPP_
#define SRC_INDEX_MAP_HPP_

// C/C++
#include <map>
#include <string>

class ParameterInput;

class IndexMap {
 protected:
  //! Protected ctor access thru static member function Instance
  IndexMap() {}

 public:
  ~IndexMap();

  static IndexMap const *GetInstance();

  static IndexMap const *InitFromAthenaInput(ParameterInput *pin);

  static void Destroy();

  size_t GetVaporId(std::string name) const {
    return vapor_index_map_.at(name);
  }

  size_t GetCloudId(std::string name) const {
    return cloud_index_map_.at(name);
  }

  size_t GetChemistryId(std::string name) const {
    return chemistry_index_map_.at(name);
  }

  size_t GetTracerId(std::string name) const {
    return tracer_index_map_.at(name);
  }

  size_t GetSpeciesId(std::string category_name) const;

 private:
  std::map<std::string, size_t> vapor_index_map_;
  std::map<std::string, size_t> cloud_index_map_;
  std::map<std::string, size_t> chemistry_index_map_;
  std::map<std::string, size_t> tracer_index_map_;
  std::map<std::string, size_t> static_index_map_;
  std::map<std::string, size_t> particle_index_map_;

  //! Pointer to the single Application instance
  static IndexMap *myindex_map_;
};

#endif  // SRC_INDEX_MAP_HPP_
