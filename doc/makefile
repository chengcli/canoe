# makefile to generate documentation
# Author: Cheng Li

SCRIPT_DIR = scripts/
EXAMPLE_DIR = examples/
DOXYGEN_DIR = doxygen/
IMG_DIR = $(DOXYGEN_DIR)html/img
MOV_DIR = $(DOXYGEN_DIR)html/mov
CF = a.tmp

EXAMPLES := straka robert
EXAMPLE_DOCS := $(addprefix $(EXAMPLE_DIR),$(EXAMPLES:%=%.dox))
PGEN_SRC := $(addprefix $(EXAMPLE_DIR),$(EXAMPLES:%=%/%.cpp))

all : start $(EXAMPLE_DOCS) doxygen
	rm -f $(CF)

start:
	@ echo 1 > $(CF)

$(EXAMPLE_DIR)%.dox: $(EXAMPLE_DIR)%
	@ mkdir -p $(EXAMPLE_DIR)
	@ perl $(SCRIPT_DIR)make_step.pl $(subst $(EXAMPLE_DIR),,$<) ../ $(shell cat $(CF)) > $@
	@ echo $(shell expr $(shell cat $(CF)) + 1) > $(CF)

doxygen: $(DOXYGEN_DIR)configure
	cd doxygen && doxygen configure
	@ mkdir -p $(IMG_DIR)
	@ mkdir -p $(MOV_DIR)
	@ ln -sf $(PWD)/images/* $(PWD)/$(IMG_DIR)
	@ ln -sf $(PWD)/movies/* $(PWD)/$(MOV_DIR)

.PHONY : start all doxygen

clean:
	rm -f $(EXAMPLE_DIR)*.dox
	rm -f $(EXAMPLE_DIR)*_plain.dox
	rm -f $(DOXYGEN_DIR)bib*.aux
	rm -f $(DOXYGEN_DIR)citelist.doc*
	rm -rf $(DOXYGEN_DIR)html
	rm -rf $(DOXYGEN_DIR)latex
