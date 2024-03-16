#! /usr/bin/env python3

# from tree_sitter import Language, Parser
from tree_sitter_languages import get_language, get_parser


def traverse(node, source_code):
    # Print the type of the node and the corresponding source code
    print(node.type)
    # print("==============")
    # print(source_code[node.start_byte:node.end_byte])
    # print("==============")

    # Recursively traverse the children of the node
    for child in node.children:
        traverse(child, source_code)


# Load the C++ language grammar
# Language.build_library(
# Store the library in the `build` directory
#  'build/my-languages.so',
# Include one or more languages
#  ['tree-sitter-cpp']
# )

# CPP_LANGUAGE = Language('build/my-languages.so', 'cpp')

# Load the language library (ensure the .so file is built and available)
# CPP_LANGUAGE = Language(tscpp.language(), 'cpp')

# Create a parser
# parser = Parser()
# parser.set_language(CPP_LANGUAGE)
parser = get_parser("cpp")

src_file = "../src/impl.hpp"

with open(src_file, "r") as file:
    code = file.read()
    tree = parser.parse(code.encode("utf-8"))

# Print the root node of the syntax tree
# print(tree.root_node.sexp())

traverse(tree.root_node, code)
