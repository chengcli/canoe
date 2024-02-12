FC = gfortran

FFLAGS = -g -std=legacy

all: crossdiskread

crossdiskread.o: crossdiskread.f
	$(FC) $(FFLAGS) -c crossdiskread.f -o $@

crossdiskread: crossdiskread.o
	$(FC) $(FFLAGS) crossdiskread.o -o $@

crossdisk-highres-combine.o: crossdisk-highres-combine.f
	$(FC) $(FFLAGS) -c crossdisk-highres-combine.f -o $@

crossdisk-highres-combine: crossdisk-highres-combine.o
	$(FC) $(FFLAGS) crossdisk-highres-combine.o -o $@

crossdisk-hirestolores.o: crossdisk-hirestolores.f
	$(FC) $(FFLAGS) -c crossdisk-hirestolores.f -o $@

crossdisk-hirestolores: crossdisk-hirestolores.o
	$(FC) $(FFLAGS) crossdisk-hirestolores.o -o $@

clean:
	/bin/rm -f *.o
