#! /usr/bin/env python3

from tree_sitter import Language, Parser

# Load the C++ language grammar
Language.build_library(
  # Store the library in the `build` directory
  'build/my-languages.so',
  # Include one or more languages
  ['tree-sitter-cpp']
)

CPP_LANGUAGE = Language('build/my-languages.so', 'cpp')

# Create a parser
parser = Parser()
parser.set_language(CPP_LANGUAGE)

src_file = "../src/impl.hpp"

with open(src_file, "r") as file:
    code = file.read()
    print(code)
    tree = parser.parse(code)

# Print the root node of the syntax tree
print(tree.root_node.sexp())
