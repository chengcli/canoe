#ifndef SRC_UTILS_FILEIO_HPP_
#define SRC_UTILS_FILEIO_HPP_

// C/C++
#include <iostream>
#include <map>
#include <string>
#include <vector>

// canoe
#include <configure.h>

// athena
#include <athena/athena.hpp>

//! test file existance
bool FileExists(std::string fname);

//! test a blank line
bool IsBlankLine(char const* line);
bool IsBlankLine(std::string const& line);

//! decomment a file
std::string DecommentFile(std::string fname);

//! get number of columns in a data table
int GetNumCols(std::string fname, char c = ' ');

//! get number of rows in a data table
int GetNumRows(std::string fname);

//! replace a character in a string
void replaceChar(char* buf, char c_old, char c_new);

template <typename T>
class AthenaArray;

char* StripLine(char* line);
char* NextLine(char* line, int num, FILE* stream);
void read_data_table(char const* fname, double** data, int* rows, int* cols);
void ReadDataTable(AthenaArray<Real>* data, std::string fname, char c = ' ');

using DataVector = std::map<std::string, std::vector<double>>;
DataVector read_data_vector(std::string fname);

#endif  // SRC_UTILS_FILEIO_HPP_
