// C/C++ header
#include <cctype>  // isspace
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

// application
#include <application/exceptions.hpp>

// utils header
#include "fileio.hpp"
#include "vectorize.hpp"

bool FileExists(std::string fname) {
  std::ifstream ifile(fname.c_str());
  return ifile.is_open();
}

bool IsBlankLine(char const* line) {
  for (char const* cp = line; *cp; ++cp) {
    if (!std::isspace(*cp)) return false;
  }
  return true;
}

bool IsBlankLine(std::string const& line) { return IsBlankLine(line.c_str()); }

std::string DecommentFile(std::string fname) {
  std::stringstream msg;
  if (!FileExists(fname)) {
    throw NotFoundError("DecommentFile", fname);
  }

  std::ifstream file(fname.c_str(), std::ios::in);
  std::string ss;
  char c;
  while (file) {
    file.get(c);
    if (c == '#') {
      while (c != '\n' && file) file.get(c);
      continue;
    }
    ss += c;
  }
  return ss;
}

int GetNumCols(std::string fname, char c) {
  std::ifstream inp(fname.c_str(), std::ios::in);
  std::string line;
  std::getline(inp, line);
  if (line.empty()) return 0;
  int cols = line[0] == c ? 0 : 1;

  for (int i = 1; i < line.length(); ++i)
    if (line[i - 1] == c && line[i] != c) cols++;
  return cols;
}

int GetNumRows(std::string fname) {
  std::ifstream inp(fname.c_str(), std::ios::in);
  std::string line;
  int rows = 0;

  while (std::getline(inp, line)) ++rows;
  return rows;
}

void replaceChar(char* buf, char c_old, char c_new) {
  int len = strlen(buf);
  for (int i = 0; i < len; ++i)
    if (buf[i] == c_old) buf[i] = c_new;
}

char* StripLine(char* line) {
  char* p = line;
  int len = strlen(line);
  // strip newline or carriage rtn
  while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
    line[--len] = 0;
  // advance to first non-whitespace
  while (isspace(*p)) p++;
  // advance to first non-whitespace
  // skip characters aftet '#'
  char* pp = p;
  while (*pp != '#' && *pp) pp++;
  *pp = 0;
  return p;
}

char* NextLine(char* line, int num, FILE* stream) {
  char* p;
  while (fgets(line, num, stream) != NULL) {
    p = StripLine(line);
    if (strlen(p) > 0) break;
  }
  return p;
}

DataVector read_data_vector(std::string fname) {
  DataVector amap;

  std::ifstream input(fname.c_str(), std::ios::in);
  std::stringstream ss, msg;
  std::string line, sbuffer;
  std::vector<std::string> field;

  if (!input.is_open()) {
    throw NotFoundError("read_data_vector", fname);
  }

  getline(input, line);
  ss.str(line);
  while (!ss.eof()) {
    ss >> sbuffer;
    field.push_back(sbuffer);
  }
  ss.clear();

  while (getline(input, line)) {
    if (line.empty()) continue;
    ss.str(line);
    for (std::vector<std::string>::iterator f = field.begin(); f != field.end();
         ++f) {
      double value;
      ss >> value;
      amap[*f].push_back(value);
    }
    ss.clear();
  }

  return amap;
}
