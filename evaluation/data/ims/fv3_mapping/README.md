(1) Code to create file with indexes used to map from ASCII IMS file onto FV3 grid. 

To compile on hera, use build.sh. 

requires following modules: 

module use /scratch1/NCEPDEV/nems/role.epic/spack-stack/spack-stack-1.6.0/envs/unified-env-rocky8/install/modulefiles/Core

module load stack-intel/2021.5.0
module load stack-intel-oneapi-mpi/2021.5.1
module load netcdf/4.7.0
module load netcdf-hdf5parallel/4.7.4

To run, submit executable (in big-mem queue)

There is a testcase and example submission script for C48 on hera at: 
/scratch2/BMC/gsienkf/Clara.Draper/DA_test_cases/snow/IMSobsproc/get_index/

(2) Use Create_nmlFile.sh to create nml file for different IMS resolutions and FV3 grids

Mike Barlage, Clara Draper, Youlong Xia.
