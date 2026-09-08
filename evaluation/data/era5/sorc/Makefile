# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .f90

 COMPILERF90    =       ifort -132
 F90FLAGS       =       -r8 -mcmodel=large
 NETCDFMOD      =       -I/apps/netcdf/4.7.0/intel/18.0.5.274/include
 NETCDFLIB      =       -L/apps/netcdf/4.7.0/intel/18.0.5.274/lib -lnetcdf -lnetcdff


OBJS =	regrid_ERA5snow_toGaussianGrid.o regrid_ERA5snow_toGaussianGrid_pt.o regrid_ERA5snow_to1Dvector.o regrid_ERA5t2m_toGaussianGrid.o

all:	regrid_ERA5snow_toGaussianGrid.exe regrid_ERA5snow_toGaussianGrid_pt.exe regrid_ERA5snow_to1Dvector.exe regrid_ERA5t2m_toGaussianGrid.exe

.f90.o:
	$(COMPILERF90) -c $(F90FLAGS) $(NETCDFMOD) $(*).f90

regrid_ERA5snow_to1Dvector.exe: $(OBJS)
	$(COMPILERF90) -o $(@) $(F90FLAGS) regrid_ERA5snow_to1Dvector.o $(NETCDFLIB)

regrid_ERA5snow_toGaussianGrid.exe: $(OBJS)
	$(COMPILERF90) -o $(@) $(F90FLAGS) regrid_ERA5snow_toGaussianGrid.o $(NETCDFLIB)

regrid_ERA5snow_toGaussianGrid_pt.exe: $(OBJS)
	$(COMPILERF90) -o $(@) $(F90FLAGS) regrid_ERA5snow_toGaussianGrid_pt.o $(NETCDFLIB)

regrid_ERA5t2m_toGaussianGrid.exe: $(OBJS)
	$(COMPILERF90) -o $(@) $(F90FLAGS) regrid_ERA5t2m_toGaussianGrid.o $(NETCDFLIB)


clean:
	rm -f *.o *.mod *.exe

#
# Dependencies:
#
