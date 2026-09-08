# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .f90

COMPILERF90    =       ifort -132 -qopenmp
F90FLAGS       =       -r8 -mcmodel=large 
NETCDFMOD      =       -I/apps/netcdf/4.7.0/intel/18.0.5.274/include
NETCDFLIB      =       -L/apps/netcdf/4.7.0/intel/18.0.5.274/lib -lnetcdf -lnetcdff

OBJS =	create_land_mapping.o

all:	create_land_mapping.exe

.f90.o:
	$(COMPILERF90) -c $(F90FLAGS) $(NETCDFMOD) $(*).f90

create_land_mapping.exe: create_land_mapping.o
	$(COMPILERF90) -o $(@) $(F90FLAGS) create_land_mapping.o $(NETCDFLIB)

clean:
	rm -f *.exe *.o

