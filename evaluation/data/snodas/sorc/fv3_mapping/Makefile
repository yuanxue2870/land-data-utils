# Makefile 
# 

.SUFFIXES:
.SUFFIXES: .o .f90

 COMPILERF90    =       ifort
 FREESOURCE     =       #-ffree-form  -ffree-line-length-none
 F90FLAGS       =       -r8 -mcmodel=large
 NETCDFMOD      =       -I/apps/netcdf/4.7.0/intel/18.0.5.274/include
 NETCDFLIB      =       -L/apps/netcdf/4.7.0/intel/18.0.5.274/lib -lnetcdf -lnetcdff

OBJS =	create_fv3_mapping.o

all:	create_fv3_mapping.exe

.f90.o:
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(NETCDFMOD) $(*).f90

create_fv3_mapping.exe: $(OBJS)
	$(COMPILERF90) -o $(@) $(F90FLAGS) $(FREESOURCE) $(NETCDFMOD) $(OBJS) $(NETCDFLIB)

clean:
	rm -f *.o *.mod *.exe

#
# Dependencies:
#
