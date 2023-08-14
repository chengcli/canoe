#ifndef SRC_OUTPUTS_OUTPUT_UTILS_HPP_
#define SRC_OUTPUTS_OUTPUT_UTILS_HPP_

// C/C++
#include <string>
#include <vector>

// athena
#include <athena/athena.hpp>

int get_num_variables(std::string grid, AthenaArray<Real> const& data);

class MetadataTable {
 protected:
  using StringTable = std::vector<std::vector<std::string>>;

  //! Protected ctor access thru static member function Instance
  MetadataTable();

 public:
  ~MetadataTable();

  static MetadataTable const* GetInstance();

  static void Destroy();

  std::string GetGridType(std::string name) const;

  std::string GetUnits(std::string name) const;

  std::string GetLongName(std::string name) const;

 private:
  StringTable table_;

  //! Pointer to the single MetadataTable instance
  static MetadataTable* myptr_;
};

#endif  // SRC_OUTPUTS_OUTPUT_UTILS_HPP_
