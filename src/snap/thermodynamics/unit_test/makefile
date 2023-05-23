#
#  Copyright (c) 2020 by Cheng Li
#

EXE = test_molecule
OBJ = $(EXE).o hydrogen.o helium.o methane.o

all : $(EXE)

$(EXE): $(OBJ)
	g++ -o $@ $(OBJ)

%.o: ../%.cpp
	g++ -c $< -o $@

clean:
	rm -f $(EXE)
	rm -f *.o

.PHONY: all clean
