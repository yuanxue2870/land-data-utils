# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .f90

COMPILERF90    =       ifort -132
F90FLAGS       =       -r8 -mcmodel=large 
#F90FLAGS       =
NETCDFMOD      =       -I${NETCDF}/include
NETCDFLIB      =       -L${NETCDF}/lib -lnetcdf -lnetcdff

OBJS =	regrid_ims_cgrid.o regrid_ims_to_gaussian.o

all:	regrid_ims_cgrid.exe regrid_ims_to_gaussian.exe

.f90.o:
	$(COMPILERF90) -c $(F90FLAGS) $(NETCDFMOD) $(*).f90

regrid_ims_cgrid.exe: regrid_ims_cgrid.o
	$(COMPILERF90) -o $(@) $(F90FLAGS) regrid_ims_cgrid.o $(NETCDFLIB)

regrid_ims_to_gaussian.exe: regrid_ims_to_gaussian.o
	$(COMPILERF90) -o $(@) $(F90FLAGS) regrid_ims_to_gaussian.o $(NETCDFLIB)

clean:
	rm -f *.exe *.o

