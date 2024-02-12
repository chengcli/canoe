#! /usr/bin/env python3
import re, sys

def stoichiometry(name):
    name_copy = name
    name = name.replace('++',' + ')
    name = re.sub(r'\+([0-9])*', r'+ \1', name)
    name = name.replace('E','E-')
    name = name.replace(' + ',' ')
    name = name.replace('=',' = ')
    name = name.split()
    eqindex = name.index('=')
    reactant, product = {}, {}
    for i in range(len(name)):
        if name[i][0].isdigit():
            num, speci = int(name[i][0]), name[i][1:]
        else:
            num, speci = 1, name[i]
        if speci == '=': continue
        if num == 0: continue
        if i < eqindex:
            reactant[speci] = num
        elif i > eqindex:
            product[speci] = num
    return reactant, product

def separate_into_chunks(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

    # Split the content into chunks based on empty lines
    chunks = [chunk.strip() for chunk in content.split('\n\n') if chunk.strip()]

    return chunks

def remove_parentheses(text):
    # Use regular expression to remove anything inside parentheses
    cleaned_text = re.sub(r'\([^)]*\)', '', text)

    return cleaned_text

def get_header_content(chunk):
    # Split the chunk into lines and get the first line
    lines = chunk.split('\n')
    header = lines[0].strip()
    content = remove_parentheses('\n'.join(lines[1:]))

    return header, content

# The first 60 characters of the header are the reaction string
def get_reaction_and_auxiliary(header):
    return header[:59].strip(), header[59:].strip()

def main(species, cross_file):
    # Separate the file content into chunks
    chunks = separate_into_chunks(cross_file)

    # Print the separated chunks
    for i, chunk in enumerate(chunks, start=1):
        header, content = get_header_content(chunk)
        reaction, auxiliary = get_reaction_and_auxiliary(header)
        reactant, product = stoichiometry(reaction)
        product_str = ' '.join([f"{speci}:{product[speci]}" for speci in product])
        if re.search(r'^%s\s*=' % species, reaction):
            #print(f" {header}\n{content}\n")
            print(f" {product_str:<60}{auxiliary}\n{content}\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python filter_cross.py species cross_file")
        print("       default cross_file is crossdisk2-use.inp-moses")
        sys.exit(1)

    species = sys.argv[1]

    if len(sys.argv) == 2:
        cross_file = "crossdisk2-use.inp-moses"
    else:
        cross_file = sys.argv[2]

    main(species, cross_file)
