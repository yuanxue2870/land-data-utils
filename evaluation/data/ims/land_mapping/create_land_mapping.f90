program create_land_mapping

use netcdf
implicit none

character*256 :: ims_regrid_filename
character*256 :: land_static_filename
character*256 :: output_filename

integer, allocatable :: data_tile      (:,:)
integer, allocatable :: data_tile_i    (:,:)
integer, allocatable :: data_tile_j    (:,:)
integer, allocatable :: model_tile     (:)
integer, allocatable :: model_tile_i   (:)
integer, allocatable :: model_tile_j   (:)
integer, allocatable :: ims_location (:,:)

integer :: error, ncid, dimid, varid, dimid_ims_i, dimid_ims_j, ims_i, ims_j, location
integer :: nlocations, idim_ims_length, jdim_ims_length, io, ierr
logical :: found, file_exists

namelist/land_mapping_nml/ ims_regrid_filename, land_static_filename, output_filename
!====================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Setup inputs and read namelist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! read namelist
 inquire(file='land_mapping.nml', exist=file_exists)

if (.not. file_exists) then
        print *, 'namelistfile does not exist, exiting'
        stop 10
endif

open (action='read', file='land_mapping.nml', iostat=ierr, newunit=io)
read (nml=land_mapping_nml, iostat=ierr, unit=io)
close (io)

! read land static file

error = nf90_open(trim(land_static_filename), NF90_NOWRITE, ncid)
  call netcdf_err(error, 'opening file: '//trim(land_static_filename) )
  
error = nf90_inq_dimid(ncid, "location", dimid)
  call netcdf_err(error, 'inquire location dimension' )

error = nf90_inquire_dimension(ncid, dimid, len = nlocations)
  call netcdf_err(error, 'reading location dimension' )
   
allocate(model_tile  (nlocations))
allocate(model_tile_i(nlocations))
allocate(model_tile_j(nlocations))

error = nf90_inq_varid(ncid, 'cube_tile', varid)
 call netcdf_err(error, 'inquire cube_tile variable' )

error = nf90_get_var(ncid, varid, model_tile)
 call netcdf_err(error, 'reading cube_tile variable' )

error = nf90_inq_varid(ncid, 'cube_i', varid)
 call netcdf_err(error, 'inquire cube_i variable' )

error = nf90_get_var(ncid, varid, model_tile_i)
 call netcdf_err(error, 'reading cube_i variable' )
    
error = nf90_inq_varid(ncid, 'cube_j', varid)
 call netcdf_err(error, 'inquire cube_j variable' )

error = nf90_get_var(ncid, varid, model_tile_j)
 call netcdf_err(error, 'reading cube_j variable' )
    
error = nf90_close(ncid)
 call netcdf_err(error, 'closing file: '//trim(land_static_filename) )
    
!====================================
! read mapping file

error = nf90_open(trim(ims_regrid_filename),nf90_nowrite, ncid)
  call netcdf_err(error, 'opening file: '//trim(ims_regrid_filename) )
    
error = nf90_inq_dimid(ncid, "idim", dimid)
  call netcdf_err(error, 'inquire idim dimension' )

error = nf90_inquire_dimension(ncid, dimid, len = idim_ims_length)
  call netcdf_err(error, 'reading idim dimension' )
   
error = nf90_inq_dimid(ncid, "jdim", dimid)
  call netcdf_err(error, 'inquire jdim dimension' )

error = nf90_inquire_dimension(ncid, dimid, len = jdim_ims_length)
  call netcdf_err(error, 'reading jdim dimension' )
   
allocate(data_tile  (idim_ims_length,jdim_ims_length))
allocate(data_tile_i(idim_ims_length,jdim_ims_length))
allocate(data_tile_j(idim_ims_length,jdim_ims_length))

error = nf90_inq_varid(ncid, 'tile', varid)
 call netcdf_err(error, 'inquire tile variable' )

error = nf90_get_var(ncid, varid, data_tile)
 call netcdf_err(error, 'reading tile variable' )
    
error = nf90_inq_varid(ncid, 'tile_i', varid)
 call netcdf_err(error, 'inquire tile_i variable' )

error = nf90_get_var(ncid, varid, data_tile_i)
 call netcdf_err(error, 'reading tile_i variable' )
    
error = nf90_inq_varid(ncid, 'tile_j', varid)
 call netcdf_err(error, 'inquire tile_j variable' )

error = nf90_get_var(ncid, varid, data_tile_j)
 call netcdf_err(error, 'reading tile_j variable' )

error = nf90_close(ncid)
 call netcdf_err(error, 'closing file: '//trim(ims_regrid_filename) )

!====================================
! create mapping to land model grid

allocate(ims_location(idim_ims_length,jdim_ims_length))

do ims_i = 1, idim_ims_length
!  print*, 'Starting: ', ims_i
do ims_j = 1, jdim_ims_length
  
  found = .false.
  do location = 1, nlocations

    if(.not.found) then
      if(model_tile  (location) == data_tile  (ims_i,ims_j) .and. &
         model_tile_i(location) == data_tile_i(ims_i,ims_j) .and. &
         model_tile_j(location) == data_tile_j(ims_i,ims_j) ) then

           ims_location(ims_i,ims_j) = location
           found = .true.

      end if
    end if

  end do

!  if(.not.found) then
!    print *, 'did not find a match'
!    print *, 'model_tile: ', model_tile(location)
!    print *, 'model_tile_i: ', model_tile_i(location)
!    print *, 'model_tile_j: ', model_tile_j(location)
!    print *, 'ims_i: ', ims_i
!    print *, 'ims_j: ', ims_j
!    print *, 'ims_tile: ', data_tile(ims_i,ims_j)
!    print *, 'ims_tile_i: ', data_tile_i(ims_i,ims_j)
!    print *, 'ims_tile_j: ', data_tile_j(ims_i,ims_j)
!    stop
!  end if

end do
end do

!====================================
! write mapping to land model grid

error = nf90_create(output_filename, ior(nf90_netcdf4,nf90_classic_model), ncid)
 call netcdf_err(error, 'creating file='//trim(output_filename) )

error = nf90_def_dim(ncid, 'idim', idim_ims_length, dimid_ims_i)
 call netcdf_err(error, 'defining idim dimension' )
 
error = nf90_def_dim(ncid, 'jdim', jdim_ims_length, dimid_ims_j)
 call netcdf_err(error, 'defining jdim dimension' )
 
error = nf90_def_var(ncid, 'mapping_location', nf90_int, (/dimid_ims_i,dimid_ims_j/), varid)
 call netcdf_err(error, 'defining mapping_location' )

error = nf90_put_att(ncid, varid, "long_name", "ims 4km data mapping in the land model vector")
 call netcdf_err(error, 'adding mapping_location long name' )

error = nf90_enddef(ncid)
 call netcdf_err(error, 'defining output_filename' )

error = nf90_put_var(ncid, varid, ims_location)
 call netcdf_err(error, 'writing ims_location variable')

error = nf90_close(ncid)
 call netcdf_err(error, 'closing output_filename')

end program
 
subroutine netcdf_err( err, string )
    
!--------------------------------------------------------------
! if a netcdf call returns an error, print out a message
! and stop processing.
!--------------------------------------------------------------
    
use netcdf

implicit none
    
integer, intent(in) :: err
character(len=*), intent(in) :: string
character(len=80) :: errmsg
    
if( err == nf90_noerr )return
errmsg = nf90_strerror(err)
print*,''
print*,'fatal error: ', trim(string), ': ', trim(errmsg)
print*,'stop.' 
stop 10

return

end subroutine netcdf_err

