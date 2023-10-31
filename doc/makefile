# makefile to generate documentation
# Author: Cheng Li

SCRIPT_DIR = scripts/
EXAMPLE_DIR = examples/
PGEN_DIR = ../examples/2019-Li-snap/
DOXYGEN_DIR = doxygen/
IMG_DIR = $(DOXYGEN_DIR)html/img
MOV_DIR = $(DOXYGEN_DIR)html/mov
CF = a.tmp

EXAMPLES := straka robert
EXAMPLE_DOCS := $(addprefix $(EXAMPLE_DIR),$(EXAMPLES:%=%.dox))
PLAIN_DOCS := $(patsubst 0d.%,%,$(EXAMPLES))
PLAIN_DOCS := $(patsubst 1d.%,%,$(PLAIN_DOCS))
PLAIN_DOCS := $(patsubst 2d.%,%,$(PLAIN_DOCS))
PLAIN_DOCS := $(patsubst 3d.%,%,$(PLAIN_DOCS))
PLAIN_DOCS := $(addsuffix _plain.txt,$(addprefix $(EXAMPLE_DIR),$(PLAIN_DOCS)))

all : start $(EXAMPLE_DOCS) $(PLAIN_DOCS) doxygen
	rm -f $(CF)

start:
	@ echo 1 > $(CF)

$(EXAMPLE_DIR)%.dox: ../$(EXAMPLE_DIR)%
	@ mkdir -p $(EXAMPLE_DIR)
	perl $(SCRIPT_DIR)make_step.pl $(subst ../$(EXAMPLE_DIR),,$<) ../ $(shell cat $(CF)) > $@
	@ echo $(shell expr $(shell cat $(CF)) + 1) > $(CF)

$(EXAMPLE_DIR)%_plain.txt: $(PGEN_DIR)%.cpp
	@ mkdir -p $(EXAMPLE_DIR)
	perl $(SCRIPT_DIR)program2plain $< > $@

doxygen: $(DOXYGEN_DIR)configure
	cd doxygen && doxygen configure
	@ mkdir -p $(IMG_DIR)
	@ mkdir -p $(MOV_DIR)
	ln -s $(PWD)/images/* $(PWD)/$(IMG_DIR)
	ln -s $(PWD)/movies/* $(PWD)/$(MOV_DIR)

.PHONY : start all doxygen

clean:
	rm -f $(EXAMPLE_DIR)*_plain.txt
	rm -f $(DOXYGEN_DIR)bib*.aux
	rm -f $(DOXYGEN_DIR)citelist.doc*
	rm -rf $(DOXYGEN_DIR)html
	rm -rf $(DOXYGEN_DIR)latex
