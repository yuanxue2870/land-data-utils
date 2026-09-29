Code uses IMS FV3-mapping index files to create land mapping files for processing 1D vector regrid output (for use of snow DA and ufs-land only run). Input is IMS FV3-mapping index file and output is land maaping file. This is openMP F90 code for 4km C1152, 1km C384, and 1km C768. it does not work for 1km C1152 as one processor cannot finish creation of land mapping suing a 8-hour limit.

To compile, do "make"
To run, submit executable (in big-mem queue) with submit_mapping.sh

Michael Barlage and Youlong Xia.
13 February 2024  
