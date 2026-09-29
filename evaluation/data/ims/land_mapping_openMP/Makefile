# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .f90

COMPILERF90    =       ifort -132 -qopenmp
F90FLAGS       =       -r8 -mcmodel=large 
#F90FLAGS       =
NETCDFMOD      =       -I${NETCDF}/include
NETCDFLIB      =       -L${NETCDF}/lib -lnetcdf -lnetcdff

OBJS =	create_land_mapping.o

all:	create_land_mapping.exe

.f90.o:
	$(COMPILERF90) -c $(F90FLAGS) $(NETCDFMOD) $(*).f90

create_land_mapping.exe: create_land_mapping.o
	$(COMPILERF90) -o $(@) $(F90FLAGS) create_land_mapping.o $(NETCDFLIB)

clean:
	rm -f *.exe *.o

