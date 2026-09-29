#!/bin/sh

# need these modules from the calling shell.
#module load intel
#module load netcdf/4.7.0 # paths hard-coded into Makefile

#cd sorc/

make -f Makefile clean 
make -f Makefile

#cd .. 

 

