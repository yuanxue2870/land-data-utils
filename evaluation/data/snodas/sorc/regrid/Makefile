# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .f90

COMPILERF90    =       ifort -132
F90FLAGS       =       -r8 -mcmodel=large 
#F90FLAGS       =
NETCDFMOD      =       -I/apps/netcdf/4.7.0/intel/18.0.5.274/include
NETCDFLIB      =       -L/apps/netcdf/4.7.0/intel/18.0.5.274/lib -lnetcdf -lnetcdff

OBJS =	regrid_snodas_cgrid.o regrid_snodas_to_gaussian.o

all:	regrid_snodas_cgrid.exe regrid_snodas_to_gaussian.exe

.f90.o:
	$(COMPILERF90) -c $(F90FLAGS) $(NETCDFMOD) $(*).f90

regrid_snodas_cgrid.exe: regrid_snodas_cgrid.o
	$(COMPILERF90) -o $(@) $(F90FLAGS) regrid_snodas_cgrid.o $(NETCDFLIB)

regrid_snodas_to_gaussian.exe: regrid_snodas_to_gaussian.o
	$(COMPILERF90) -o $(@) $(F90FLAGS) regrid_snodas_to_gaussian.o $(NETCDFLIB)

clean:
	rm -f *.exe *.o

