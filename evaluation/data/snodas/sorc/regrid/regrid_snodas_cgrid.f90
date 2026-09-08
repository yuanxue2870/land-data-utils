program regrid_snodas_cgrid

use netcdf

implicit none

character*256 :: land_static_filename
character*256 :: snodas_regrid_filename
character*256 :: snodas_filename 
character*256 :: output_filename 

character*256 :: snodas_filepath
character*256 :: regrid_outpath

character*10 :: obs_source, fv3_grid

double precision  :: sec_since
integer           :: offset_ss
integer           :: dim_time, id_time, id_lat, id_lon
character*19      :: current_date  ! current date
character*19      :: since_date = "1970-01-01 00:00:00"
character*8       :: datestring
character*4       :: yyyy
character*2       :: mm, dd

real     , allocatable :: liquid_precip(:,:)
real     , allocatable :: solid_precip(:,:)
real     , allocatable :: snowdepth(:,:)
real     , allocatable :: swe(:,:)

integer  , allocatable :: mapping_location(:,:)

real     , allocatable :: latitude(:)
real     , allocatable :: longitude(:)

real     , allocatable :: data_regrid        (:)

logical :: file_exists

integer :: error, ncid, dimid, varid(4), io, ierr
integer :: nlocations, idim_snodas_length, jdim_snodas_length

real*4  fillVal;
fillVal = -9999.0

namelist/snodas_regrid_nml/ obs_source, fv3_grid, land_static_filename, snodas_regrid_filename, snodas_filepath, regrid_outpath 

! read namelist

inquire(file='snodas_regrid.nml', exist=file_exists)

if (.not. file_exists) then
    print *, 'namelistfile does not exist, exiting'
    stop 10
endif

open (action='read', file='snodas_regrid.nml', iostat=ierr, newunit=io)
read (nml=snodas_regrid_nml, iostat=ierr, unit=io)
close (io)

open(18, file='date_input.txt', status='old')
read(18,'(A4)') yyyy
read(18,'(A2)') mm
read(18,'(A2)') dd
close(18)

offset_ss=0
datestring=yyyy//mm//dd
current_date=yyyy//'-'//mm//'-'//dd//' 18:00:00'
call calc_sec_since(since_date, current_date, offset_ss, sec_since)

!====================================
! read land static file

error = nf90_open(trim(land_static_filename), NF90_NOWRITE, ncid)
  call netcdf_err(error, 'opening file: '//trim(land_static_filename) )
  
error = nf90_inq_dimid(ncid, "location", dimid)
  call netcdf_err(error, 'inquire location dimension' )

error = nf90_inquire_dimension(ncid, dimid, len = nlocations)
  call netcdf_err(error, 'reading location dimension' )

allocate(latitude(nlocations))
allocate(longitude(nlocations))   

error = nf90_inq_varid(ncid, 'latitude', varid(1))
 call netcdf_err(error, 'inquire latitude variable' )
error = nf90_get_var(ncid, varid(1), latitude)
 call netcdf_err(error, 'reading latitude variable' )

error = nf90_inq_varid(ncid, 'longitude', varid(1))
 call netcdf_err(error, 'inquire longitude variable' )
error = nf90_get_var(ncid, varid(1), longitude)
 call netcdf_err(error, 'reading longitude variable' )

error = nf90_close(ncid)
 call netcdf_err(error, 'closing file: '//trim(land_static_filename) )

!====================================
! read mapping file

error = nf90_open(trim(snodas_regrid_filename),nf90_nowrite, ncid)
  call netcdf_err(error, 'opening file: '//trim(snodas_regrid_filename) )
    
error = nf90_inq_dimid(ncid, "idim", dimid)
  call netcdf_err(error, 'inquire idim dimension' )

error = nf90_inquire_dimension(ncid, dimid, len = idim_snodas_length)
  call netcdf_err(error, 'reading idim dimension' )
   
error = nf90_inq_dimid(ncid, "jdim", dimid)
  call netcdf_err(error, 'inquire jdim dimension' )

error = nf90_inquire_dimension(ncid, dimid, len = jdim_snodas_length)
  call netcdf_err(error, 'reading jdim dimension' )

allocate(mapping_location(idim_snodas_length,jdim_snodas_length))

error = nf90_inq_varid(ncid, 'mapping_location', varid(1))
 call netcdf_err(error, 'inquire mapping_location variable' )

error = nf90_get_var(ncid, varid(1), mapping_location)
 call netcdf_err(error, 'reading mapping_location variable' )
    
error = nf90_close(ncid)
 call netcdf_err(error, 'closing file: '//trim(snodas_regrid_filename) )

!====================================
! read data file
allocate(snowdepth(idim_snodas_length,jdim_snodas_length))
allocate(swe(idim_snodas_length,jdim_snodas_length))
allocate(liquid_precip(idim_snodas_length,jdim_snodas_length))
allocate(solid_precip(idim_snodas_length,jdim_snodas_length))

snodas_filename = trim(snodas_filepath)//trim(obs_source)//'_unmasked_'//trim(datestring)//'.nc'

!write(*,*) 'snodas_filename: ', snodas_filename

error = nf90_open(trim(snodas_filename),nf90_nowrite, ncid)
  call netcdf_err(error, 'opening file: '//trim(snodas_filename) )

error = nf90_inq_varid(ncid, 'liquid_precipitation', varid(1))
 call netcdf_err(error, 'inquire liquid_precipitation variable' )

error = nf90_get_var(ncid, varid(1), liquid_precip)
 call netcdf_err(error, 'reading liquid_precipitation variable' )

error = nf90_inq_varid(ncid, 'solid_precipitation', varid(1))
 call netcdf_err(error, 'inquire solid_precipitation variable' )

error = nf90_get_var(ncid, varid(1), solid_precip)
 call netcdf_err(error, 'reading solid_precip variable' )
   
error = nf90_inq_varid(ncid, 'snow_depth', varid(1))
 call netcdf_err(error, 'inquire snow_depth variable' )

error = nf90_get_var(ncid, varid(1), snowdepth)
 call netcdf_err(error, 'reading snow_depth variable' )
    
error = nf90_inq_varid(ncid, 'snow_water_equivalent', varid(1))
 call netcdf_err(error, 'inquire snow_water_equivalent variable' )

error = nf90_get_var(ncid, varid(1), swe)
 call netcdf_err(error, 'reading snow_water_equivalent variable' )

error = nf90_close(ncid)
 call netcdf_err(error, 'closing file: '//trim(snodas_filename) )

allocate(data_regrid(nlocations))

!====================================
! write remapped fields

output_filename = trim(regrid_outpath)//trim(fv3_grid)//'/'//trim(obs_source)//'_snow.'//trim(fv3_grid)//'.'//trim(datestring)//'.nc'

!write(*,*) output_filename

error = nf90_create(output_filename, ior(nf90_netcdf4,nf90_classic_model), ncid)
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
error = nf90_put_att(ncid, id_lon, "units", "degrees_east")
  call netcdf_err(error, 'defining lon units' )

!--- define latitude
error = nf90_def_var(ncid, 'lat', nf90_double, dimid, id_lat)
  call netcdf_err(error, 'defining lat' )
error = nf90_put_att(ncid, id_lat, "long_name", "latitude")
  call netcdf_err(error, 'defining lat long name' )
error = nf90_put_att(ncid, id_lat, "units", "degrees_north")
  call netcdf_err(error, 'defining lat units' )

!======================================
! liquid precipitation

error = nf90_def_var(ncid, 'liquid_precipitation', nf90_float, (/dimid, dim_time/), varid(1))
 call netcdf_err(error, 'defining liquid precipitation' )
error = nf90_put_att(ncid, varid(1), "long_name", "SNODAS liquid precipitation")
  call netcdf_err(error, 'defining liquid precipitation long name' )
error = nf90_put_att(ncid, varid(1), "units", "mm/day")
  call netcdf_err(error, 'defining liquid precipitation units' )
error = nf90_put_att(ncid, varid(1), "_FillValue", fillVal)
 call netcdf_err(error, 'adding liquid precipitation _FillValue' )

! ======= solid precipitation =============

error = nf90_def_var(ncid, 'solid_precipitation', nf90_float, (/dimid, dim_time/), varid(2))
 call netcdf_err(error, 'defining solid precipitation' )
error = nf90_put_att(ncid, varid(2), "long_name", "SNODAS solid precipitation")
  call netcdf_err(error, 'defining solid precipitation long name' )
error = nf90_put_att(ncid, varid(2), "units", "mm/day")
  call netcdf_err(error, 'defining solid precipitation units' )
error = nf90_put_att(ncid, varid(2), "_FillValue", fillVal)
 call netcdf_err(error, 'adding solid precipitation _FillValue' )

!======= snow depth ======================

error = nf90_def_var(ncid, 'snow_depth', nf90_float, (/dimid, dim_time/), varid(3))
 call netcdf_err(error, 'defining snow depth' )
error = nf90_put_att(ncid, varid(3), "long_name", "SNODAS snow depth")
  call netcdf_err(error, 'defining snow depth long name' )
error = nf90_put_att(ncid, varid(3), "units", "mm")
  call netcdf_err(error, 'defining snow depth units' )
error = nf90_put_att(ncid, varid(3), "_FillValue", fillVal)
 call netcdf_err(error, 'adding snow depth _FillValue' )

!============== snow water equivalent ========
error = nf90_def_var(ncid, 'snow_water_equivalent', nf90_float, (/dimid, dim_time/), varid(4))
 call netcdf_err(error, 'defining snow water equivalent' )
error = nf90_put_att(ncid, varid(4), "long_name", "SNODAS snow water equivalent")
  call netcdf_err(error, 'defining SWE long name' )
error = nf90_put_att(ncid, varid(4), "units", "mm")
  call netcdf_err(error, 'defining snow water equivalent units' )
 error = nf90_put_att(ncid, varid(4), "_FillValue", fillVal)
 call netcdf_err(error, 'adding snow water equivalent _FillValue' )

error = nf90_enddef(ncid)
 call netcdf_err(error, 'defining output_filename' )

! --- put time, lat, lon data
error = nf90_put_var(ncid, id_time, sec_since)
    call netcdf_err(error, 'writing time record')

error = nf90_put_var(ncid, id_lon, longitude)
    call netcdf_err(error, 'writing lon record')

error = nf90_put_var(ncid, id_lat, latitude)
  call netcdf_err(error, 'writing lat record')

call  average_snodasSnow(idim_snodas_length, jdim_snodas_length, nlocations, mapping_location, liquid_precip, data_regrid, 1)
error = nf90_put_var(ncid, varid(1), data_regrid, start = (/1,1/), count = (/nlocations,1/))
 call netcdf_err(error, 'writing snodas liquid precipitation regrid variable')

call average_snodasSnow(idim_snodas_length, jdim_snodas_length, nlocations, mapping_location, solid_precip, data_regrid, 1)
error = nf90_put_var(ncid, varid(2), data_regrid, start = (/1,1/), count = (/nlocations,1/))
 call netcdf_err(error, 'writing snodas solid precipitation regrid variable')

call  average_snodasSnow(idim_snodas_length, jdim_snodas_length, nlocations, mapping_location, snowdepth, data_regrid, 2)
error = nf90_put_var(ncid, varid(3), data_regrid, start = (/1,1/), count = (/nlocations,1/))
 call netcdf_err(error, 'writing snodas snow depth regrid variable')

call  average_snodasSnow(idim_snodas_length, jdim_snodas_length, nlocations, mapping_location, swe, data_regrid, 2)
error = nf90_put_var(ncid, varid(4), data_regrid, start = (/1,1/), count = (/nlocations,1/))
 call netcdf_err(error, 'writing snodas snow water equivalent regrid variable')

error = nf90_close(ncid)
 call netcdf_err(error, 'closing output_filename')

end program


subroutine average_snodasSnow(idim_snodas_length, jdim_snodas_length, nlocations, mapping_location, snowdepth, data_regrid, mmcal)

!--------------------------------------------------------------
! computing averaged snodas data values  for a given FV3 Cxx grid
!-------------------------------------------------------------

implicit none

integer :: nlocations, idim_snodas_length, jdim_snodas_length

real, dimension(idim_snodas_length,jdim_snodas_length) :: snowdepth
real, dimension(nlocations)                            :: data_regrid

integer  , dimension(idim_snodas_length,jdim_snodas_length) :: mapping_location
integer  , dimension(nlocations)                            :: land_grids
real     , dimension(nlocations)                            :: snodas_regrid  

integer ::  snodas_i, snodas_j, location
integer :: land_check
integer :: mmcal    ! mmcal=1, multiply scale factor 0.1; mmcal=2, keep mm unit

real*4  fillVal;
fillVal = -9999.0

snodas_regrid = 0.0
land_grids = 0

do snodas_i = 1, idim_snodas_length
do snodas_j = 1, jdim_snodas_length
  
  location = mapping_location(snodas_i,snodas_j)

  if(location > 0) then

    land_check = 0 
    
! confirm removal of missing data
 
    if(snowdepth(snodas_i,snodas_j) >=0 ) land_check = 1

    if(land_check == 1) then

      land_grids(location) = land_grids(location) + 1
      snodas_regrid(location)  =  snodas_regrid(location) +  snowdepth(snodas_i,snodas_j)

    end if ! land_check == 1

  end if ! location > 0

end do
end do

if (mmcal == 1) then
where(land_grids >  0) data_regrid = 0.1*snodas_regrid/real(land_grids)
end if

if (mmcal == 2) then
where(land_grids >  0) data_regrid = snodas_regrid/real(land_grids)
end if

where(land_grids == 0) data_regrid = fillVal

return

end subroutine average_snodasSnow

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
 
