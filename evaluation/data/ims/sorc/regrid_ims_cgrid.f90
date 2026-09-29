program regrid_ims_cgrid

use netcdf

implicit none

character*100 :: regrid_filepath
character*100 :: ims_filepath
character*100 :: regrid_outpath

character*256 :: ims_source_path
character*256 :: ims_out_path

character*256 :: ims_regrid_filename
character*256 :: land_static_filename
character*256 :: ims_filename 
character*256 :: output_filename 

character*10   :: obs_source
character*10   :: fv3_grid
character*3    :: ims_res

character*7   :: datestring
double precision  :: sec_since
integer           :: offset_ss
integer           :: dim_time, id_time, id_lat, id_lon
character*19      :: current_date  ! current date
character*19      :: since_date = "1970-01-01 00:00:00"
character*8       :: daystring
character*4       :: yyyy
character*2       :: mm, dd

integer*2, allocatable :: imssnow_flag    (:,:)
integer  , allocatable :: mapping_location(:,:)

real     , allocatable :: latitude(:)
real     , allocatable :: longitude(:)

integer*2, allocatable :: data_regrid        (:)

integer :: error, ncid, dimid, varid, io, ierr
integer :: nlocations, idim_ims_length, jdim_ims_length
integer*2, parameter ::  nodata_int = -32767

logical :: file_exists

namelist/ims_regrid_nml/ obs_source, fv3_grid, land_static_filename, regrid_filepath, ims_filepath, regrid_outpath 

! read namelist

inquire(file='ims_regrid.nml', exist=file_exists)

if (.not. file_exists) then
    print *, 'namelistfile does not exist, exiting'
    stop 10
endif

open (action='read', file='ims_regrid.nml', iostat=ierr, newunit=io)
read (nml=ims_regrid_nml, iostat=ierr, unit=io)
close (io)

ims_res=obs_source(4:6)

! create 1km and 4km ims_regrid_filename, actual ims_source and ims_out path
if (trim(obs_source)=="IMS1km") then
    ims_source_path=trim(ims_filepath)//trim(ims_res)//'/'
    ims_regrid_filename=trim(regrid_filepath)//trim(obs_source)//'_to_land_mapping.'//trim(fv3_grid)//'.nc'
    ims_out_path=trim(regrid_outpath)
elseif (trim(obs_source) == "IMS4km" ) then
    ims_source_path=trim(ims_filepath)
    ims_regrid_filename=trim(regrid_filepath)//trim(obs_source)//'_to_land_mapping.'//trim(fv3_grid)//'.nc'
    ims_out_path=trim(regrid_outpath)    
else
   write(6,*) 'obs_source not recognised', obs_source
   stop 10
endif

! read year, month, day and date sring

open(18, file='date_input.txt', status='old')
read(18,'(A7)') datestring
read(18,'(A8)') daystring
read(18,'(A4)') yyyy
read(18,'(A2)') mm
read(18,'(A2)') dd
close(18)

offset_ss=0
current_date=yyyy//'-'//mm//'-'//dd//' 18:00:00'
call calc_sec_since(since_date, current_date, offset_ss, sec_since)

!====================================
! read land static file
!write(*,*) trim(land_static_filename)

error = nf90_open(trim(land_static_filename), NF90_NOWRITE, ncid)
  call netcdf_err(error, 'opening file: '//trim(land_static_filename) )
  
error = nf90_inq_dimid(ncid, "location", dimid)
  call netcdf_err(error, 'inquire location dimension' )

error = nf90_inquire_dimension(ncid, dimid, len = nlocations)
  call netcdf_err(error, 'reading location dimension' )

allocate(latitude(nlocations))
allocate(longitude(nlocations))   

error = nf90_inq_varid(ncid, 'latitude', varid)
 call netcdf_err(error, 'inquire latitude variable' )
error = nf90_get_var(ncid, varid, latitude)
 call netcdf_err(error, 'reading latitude variable' )

error = nf90_inq_varid(ncid, 'longitude', varid)
 call netcdf_err(error, 'inquire longitude variable' )
error = nf90_get_var(ncid, varid, longitude)
 call netcdf_err(error, 'reading longitude variable' )

error = nf90_close(ncid)
 call netcdf_err(error, 'closing file: '//trim(land_static_filename) )

!====================================
! read mapping file
!write(*,*) trim(ims_regrid_filename)

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
   
allocate(mapping_location(idim_ims_length,jdim_ims_length))

error = nf90_inq_varid(ncid, 'mapping_location', varid)
 call netcdf_err(error, 'inquire mapping_location variable' )

error = nf90_get_var(ncid, varid, mapping_location)
 call netcdf_err(error, 'reading mapping_location variable' )
    
error = nf90_close(ncid)
 call netcdf_err(error, 'closing file: '//trim(ims_regrid_filename) )

!====================================
! read data file
allocate(imssnow_flag(idim_ims_length,jdim_ims_length))

ims_filename = trim(ims_source_path)//'/ims'//trim(datestring)//'_'//trim(ims_res)//'_v1.3.nc'

!write(*,*) trim(ims_filename)

error = nf90_open(trim(ims_filename),nf90_nowrite, ncid)
  call netcdf_err(error, 'opening file: '//trim(ims_filename) )
    
error = nf90_inq_varid(ncid, 'IMS_Surface_Values', varid)
 call netcdf_err(error, 'inquire IMS_Surface_Values variable' )

error = nf90_get_var(ncid, varid, imssnow_flag)
 call netcdf_err(error, 'reading imssnow_flag variable' )
    
error = nf90_close(ncid)
 call netcdf_err(error, 'closing file: '//trim(ims_filename) )

where(imssnow_flag == 0 ) imssnow_flag = nodata_int ! set outside range to NA
where(imssnow_flag == 1 ) imssnow_flag = nodata_int ! set sea to NA
where(imssnow_flag == 3 ) imssnow_flag = nodata_int ! set sea ice to NA
where(imssnow_flag == 2 ) imssnow_flag = 0          ! set land, no snow to 0
where(imssnow_flag == 4 ) imssnow_flag = 1          ! set snow on land to 1

allocate(data_regrid(nlocations))

!====================================
! write remapped fields

output_filename =trim(ims_out_path)//'/IMSscf.'//trim(fv3_grid)//'.'//ims_res//'.'//trim(daystring)//'.nc'

error = nf90_create(trim(output_filename), ior(nf90_netcdf4,nf90_classic_model), ncid)
 call netcdf_err(error, 'creating file='//trim(output_filename) )

! --- define location and time

error = nf90_def_dim(ncid, 'location', nlocations, dimid)
 call netcdf_err(error, 'defining location dimension' )
error = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dim_time)
  call netcdf_err(error, 'defining time dimension' )

error = nf90_def_var(ncid, 'time', nf90_double, dim_time, id_time)
  call netcdf_err(error, 'defining time' )
error = nf90_put_att(ncid, id_time, "long_name", "time")
  call netcdf_err(error, 'defining time long name' )
error = nf90_put_att(ncid, id_time, "units", "seconds since "//since_date)
  call netcdf_err(error, 'defining time units' )

!--- define longitude
error = nf90_def_var(ncid, 'lon', nf90_double, dimid, id_lon)
  call netcdf_err(error, 'defining lon' )
error = nf90_put_att(ncid, id_lon, "long_name", "longitude")
  call netcdf_err(error, 'defining lon long name' )

!--- define latitude
error = nf90_def_var(ncid, 'lat', nf90_double, dimid, id_lat)
  call netcdf_err(error, 'defining lat' )
error = nf90_put_att(ncid, id_lat, "long_name", "latitude")
  call netcdf_err(error, 'defining lat long name' )

!======================================
! percent of snow

error = nf90_def_var(ncid, 'snow_cover_fraction', nf90_short, (/dimid, dim_time/), varid)
 call netcdf_err(error, 'defining snow_cover_fraction' )

error = nf90_put_att(ncid, varid, "_FillValue", nodata_int)
 call netcdf_err(error, 'adding snow_cover_fraction _FillValue' )

error = nf90_enddef(ncid)
 call netcdf_err(error, 'defining output_filename' )

! --- put time, lat, lon data
error = nf90_put_var(ncid, id_time, sec_since)
    call netcdf_err(error, 'writing time record')

error = nf90_put_var(ncid, id_lon, longitude)
    call netcdf_err(error, 'writing lon record')

error = nf90_put_var(ncid, id_lat, latitude)
  call netcdf_err(error, 'writing lat record')

call  average_imsSnow(idim_ims_length, jdim_ims_length, nlocations, mapping_location, imssnow_flag, data_regrid)
error = nf90_put_var(ncid, varid, data_regrid, start = (/1,1/), count = (/nlocations,1/))
 call netcdf_err(error, 'writing imssnow_regrid variable')

error = nf90_close(ncid)
 call netcdf_err(error, 'closing output_filename')

end program


subroutine average_imsSnow(idim_ims_length, jdim_ims_length, nlocations, mapping_location, imssnow_flag, data_regrid)

!--------------------------------------------------------------
! computing averaged ims snow cover values  for a given FV3 Cxx grid
!-------------------------------------------------------------

implicit none

integer :: nlocations, idim_ims_length, jdim_ims_length

integer*2, dimension(idim_ims_length,jdim_ims_length) :: imssnow_flag
integer*2, dimension(nlocations)                      :: data_regrid

integer  , dimension(idim_ims_length,jdim_ims_length) :: mapping_location
integer  , dimension(nlocations)                          :: land_grids
real     , dimension(nlocations)                          :: ims_regrid  

integer ::  ims_i, ims_j, location
integer :: land_check
integer*2, parameter ::  nodata_int = -32767

ims_regrid = 0.0
land_grids = 0

do ims_i = 1, idim_ims_length
do ims_j = 1, jdim_ims_length
  
  location = mapping_location(ims_i,ims_j)

  if(location > 0) then

    land_check = 0 
    
! confirm removal of missing data
 
    if(imssnow_flag(ims_i,ims_j) >=0 ) land_check = 1

    if(land_check == 1) then

      land_grids(location) = land_grids(location) + 1
      ims_regrid(location)  =  ims_regrid(location) +  imssnow_flag(ims_i,ims_j)

    end if ! land_check == 1

  end if ! location > 0

end do
end do

where(land_grids >  0) ims_regrid = 100.0*real(ims_regrid)/real(land_grids)

data_regrid = nint(ims_regrid)
where(land_grids == 0) data_regrid = nodata_int

return

end subroutine average_imsSnow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate time in seconds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_sec_since(since_date, current_date, offset_ss, sec_since)

! calculate number of seconds between since_date and current_date

double precision      :: sec_since
character*19 :: since_date, current_date  ! format: yyyy-mm-dd hh:nn:ss
integer      :: offset_ss
integer      :: since_yyyy, since_mm, since_dd, since_hh, since_nn, since_ss
integer      :: current_yyyy, current_mm, current_dd, current_hh, current_nn, current_ss
logical      :: leap_year = .false.
integer      :: iyyyy, imm
integer, dimension(12), parameter :: days_in_mm = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  sec_since = 0

  read(since_date( 1: 4),  '(i4)') since_yyyy
  read(since_date( 6: 7),  '(i2)') since_mm
  read(since_date( 9:10),  '(i2)') since_dd
  read(since_date(12:13),  '(i2)') since_hh
  read(since_date(15:16),  '(i2)') since_nn
  read(since_date(18:19),  '(i2)') since_ss

  read(current_date( 1: 4),  '(i4)') current_yyyy
  read(current_date( 6: 7),  '(i2)') current_mm
  read(current_date( 9:10),  '(i2)') current_dd
  read(current_date(12:13),  '(i2)') current_hh
  read(current_date(15:16),  '(i2)') current_nn
  read(current_date(18:19),  '(i2)') current_ss

! not worrying about the complexity of non-recent leap years
! calculate number of seconds in all years
  do iyyyy = since_yyyy, current_yyyy
    if(mod(iyyyy,4) == 0) then
      sec_since = sec_since + 366*86400
    else
      sec_since = sec_since + 365*86400
    end if
  end do

! remove seconds from since_year
  if(mod(since_yyyy,4) == 0) leap_year = .true.

  do imm = 1,since_mm-1
    sec_since = sec_since - days_in_mm(imm)*86400
  end do

  if(leap_year .and. since_mm > 2) sec_since = sec_since - 86400

  sec_since = sec_since - (since_dd - 1) * 86400

  sec_since = sec_since - (since_hh) * 3600

  sec_since = sec_since - (since_nn) * 60

  sec_since = sec_since - (since_ss)

! remove seconds in current_year

  leap_year = .false.
  if(mod(current_yyyy,4) == 0) leap_year = .true.

  do imm = current_mm+1, 12
    sec_since = sec_since - days_in_mm(imm)*86400
  end do
  if(leap_year .and. current_mm < 3) sec_since = sec_since - 86400

  sec_since = sec_since - (days_in_mm(current_mm) - current_dd) * 86400

  sec_since = sec_since - (23 - current_hh) * 3600

  sec_since = sec_since - (59 - current_nn) * 60

  sec_since = sec_since - (60 - current_ss)

  sec_since = sec_since + offset_ss

end subroutine calc_sec_since

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
 
