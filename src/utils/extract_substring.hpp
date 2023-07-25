#ifndef SRC_UTILS_EXTRACT_SUBSTRING_HPP_
#define SRC_UTILS_EXTRACT_SUBSTRING_HPP_

#include <string>

// extract name before delimiter
std::string extract_first(std::string first_second, std::string delimiter);

// extract name after delimiter
std::string extract_second(std::string first_second, std::string delimiter);

#endif
