// C/C++ header
#include <cctype>  // isspace
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

// harp header
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
    msg << "### FATAL ERROR in DecommentFile. File " << fname
        << " doesn't exist." << std::endl;
    throw std::runtime_error(msg.str().c_str());
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

template <>
std::vector<std::string> Vectorize(const char* cstr, const char* delimiter) {
  std::vector<std::string> arr;
  char str[1028], *p;
  strcpy(str, cstr);
  p = std::strtok(str, delimiter);
  while (p != NULL) {
    arr.push_back(p);
    p = std::strtok(NULL, delimiter);
  }
  return arr;
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
