// C/C++
#include <string>

// extract name before delimiter
std::string extract_first(std::string first_second, std::string delimiter) {
  size_t delimiter_pos = first_second.find(delimiter);
  return first_second.substr(0, delimiter_pos);
}

// extract name after delimiter
std::string extract_second(std::string first_second, std::string delimiter) {
  size_t delimiter_pos = first_second.find(delimiter);
  return first_second.substr(delimiter_pos + delimiter.length());
}
