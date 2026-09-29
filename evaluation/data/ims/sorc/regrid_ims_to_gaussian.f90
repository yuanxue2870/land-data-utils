program regrid_ims_gaussian

use netcdf

implicit none

character*100 :: map_filepath
character*100 :: ims_filepath
character*100 :: regrid_outpath

character*256 :: ims_source_path
character*256 :: ims_out_path

character*256 :: land_fix_filename
character*256 :: ims_filename
character*256 :: map_filename 
character*256 :: output_filename 

character*10   :: obs_source
character*10   :: fv3_grid
character*3    :: ims_res

character*7   :: datestring
double precision  :: sec_since
integer           :: offset_ss
integer           :: dim_time, id_time, id_lat, id_lon
integer           :: dim_id_j, dim_id_i
integer           :: i, j, itotal, tile_i, tile_j 
character*19      :: current_date  ! current date
character*19      :: since_date = "1970-01-01 00:00:00"
character*8       :: daystring
character*4       :: yyyy
character*2       :: mm, dd

integer*2, allocatable :: imssnow_flag    (:,:)
integer*2, allocatable :: lookup_i        (:,:)
integer*2, allocatable :: lookup_j        (:,:)

real     , allocatable :: latitude(:)
real     , allocatable :: longitude(:)
real     , allocatable :: xcb(:)  ! destination lon
real     , allocatable :: ycb (:) ! destination lat
real     , allocatable :: lon2d(:,:)
real     , allocatable :: lat2d(:,:)

real     , allocatable:: pt_snow_cover(:,:)
real     , allocatable:: pt_land(:,:)
real     , allocatable:: pt_snow(:,:)

integer :: error, ncid, dimid, varid, io, ierr
integer :: source_i_size,source_j_size
integer :: destination_lons, destination_lats, destination_locs 
integer*2, parameter ::  nodata_int = -32767
logical :: file_exists

namelist/ims_regrid_nml/ destination_lons, destination_lats, obs_source, fv3_grid, land_fix_filename,  map_filepath, ims_filepath, regrid_outpath 

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
destination_locs=destination_lats*destination_lons

allocate(xcb(destination_locs))
allocate(ycb(destination_locs))
allocate(latitude(destination_lats))
allocate(longitude(destination_lons))
allocate(lon2d(destination_lons,destination_lats))
allocate(lat2d(destination_lons,destination_lats))

allocate(pt_snow_cover(destination_lons,destination_lats))
allocate(pt_land(destination_lons,destination_lats))
allocate(pt_snow(destination_lons,destination_lats))

! create 1km and 4km ims_regrid_filename, actual ims_source and ims_out path
if (trim(obs_source)=="IMS1km") then
    source_i_size=24576
    source_j_size=24576
    ims_source_path=trim(ims_filepath)//trim(ims_res)//'/'
    ims_out_path=trim(regrid_outpath)//trim(obs_source)//'/'//trim(fv3_grid)
elseif (trim(obs_source) == "IMS4km" ) then
    source_i_size=6144
    source_j_size=6144
    ims_source_path=trim(ims_filepath)
    ims_out_path=trim(regrid_outpath)//trim(obs_source)//'/'//trim(fv3_grid)    
else
   write(6,*) 'obs_source not recognised', obs_source
   stop 10
endif

allocate(imssnow_flag(source_i_size,source_j_size))
allocate(lookup_i(source_i_size,source_j_size))
allocate(lookup_j(source_i_size,source_j_size))

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
! read land fix gassian file
!write(*,*) trim(land_fix_filename)

error = nf90_open(trim(land_fix_filename), NF90_NOWRITE, ncid)
  call netcdf_err(error, 'opening file: '//trim(land_fix_filename) )
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read lattitude and longitude data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

error = nf90_inq_varid(ncid, "yc_b", varid)
 call netcdf_err(error, 'inquire latitude variable' )
error = nf90_get_var(ncid, varid, ycb)
 call netcdf_err(error, 'reading latitude variable' )

error =  nf90_inq_varid(ncid, "xc_b", varid)
 call netcdf_err(error, 'inquire longitude variable' )
error = nf90_get_var(ncid, varid, xcb)
 call netcdf_err(error, 'reading longitude variable' )

error = nf90_close(ncid)
 call netcdf_err(error, 'closing file: '//trim(land_fix_filename) )

!==============================================================
! Convert 1D vector into 2D array for soil moisture and error
! Keep lats and lons as 1D to save storage space
!===============================================================
 itotal=0
 do j = 1, destination_lats
 do i = 1, destination_lons
   itotal=itotal+1
   lat2d(i,j) = ycb(itotal)
   lon2d(i,j) = xcb(itotal)
 end do
 end do

 do j = 1,destination_lats
    latitude(j) = lat2d(1,j)
 end do

 do i = 1, destination_lons
   longitude(i) = lon2d(i,1)
 end do
!====================================
! read mapping file

map_filename=trim(map_filepath)//trim(obs_source)//'_'//trim(fv3_grid)//'_mapping.nc'
!write(*,*) trim(map_filename)

error = nf90_open(trim(map_filename),nf90_nowrite, ncid)
  call netcdf_err(error, 'opening file: '//trim(map_filename) )

error = nf90_inq_varid(ncid, "pt_i_ims", varid)
  call netcdf_err(error, 'inquire pt_i_ims variable' ) 
error = nf90_get_var(ncid, varid , lookup_i)
  call netcdf_err(error, 'reading pt_i_ims variable' ) 

error = nf90_inq_varid(ncid, "pt_j_ims", varid)
  call netcdf_err(error, 'inquire pt_j_ims variable' )
error = nf90_get_var(ncid, varid , lookup_j)
  call netcdf_err(error, 'reading pt_j_ims variable' )

error = nf90_close(ncid)
 call netcdf_err(error, 'closing file: '//trim(map_filename) )

!====================================
! read IMS NetCDF data file

ims_filename = trim(ims_source_path)//yyyy//'/ims'//trim(datestring)//'_'//trim(ims_res)//'_v1.3.nc'

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

! calculate fraction

pt_snow         = 0
pt_land         = 0
pt_snow_cover   = -9999.0

do j=1, source_j_size
 do i=1, source_i_size

  if(imssnow_flag(i,j) >= 0) then
     tile_i = lookup_i(i,j)
     tile_j = lookup_j(i,j)
     pt_land(tile_i,tile_j)  = pt_land(tile_i,tile_j) + 1 
 
     if(imssnow_flag(i,j)==1) then
        pt_snow(tile_i,tile_j) = pt_snow(tile_i,tile_j) + 1
     end if
   end if
 end do
enddo

where(pt_land > 0) pt_snow_cover = 100.0*(pt_snow/pt_land)

!====================================
! write remapped fields

output_filename =trim(ims_out_path)//'/IMSscf.'//trim(fv3_grid)//'.'//ims_res//'.'//trim(daystring)//'.nc'

error = nf90_create(trim(output_filename), ior(nf90_netcdf4,nf90_classic_model), ncid)
 call netcdf_err(error, 'creating file='//trim(output_filename) )

! --- define dimensions and time

error = nf90_def_dim(ncid, "lat"   , destination_lats , dim_id_j)
 call netcdf_err(error, 'defining latitude dimension' )

error =  nf90_def_dim(ncid, "lon"   , destination_lons , dim_id_i)
 call netcdf_err(error, 'defining longitude dimension' )

error = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dim_time)
  call netcdf_err(error, 'defining time dimension' )

error = nf90_def_var(ncid, 'time', nf90_double, dim_time, id_time)
  call netcdf_err(error, 'defining time' )
error = nf90_put_att(ncid, id_time, "long_name", "time")
  call netcdf_err(error, 'defining time long name' )
error = nf90_put_att(ncid, id_time, "units", "seconds since "//since_date)
  call netcdf_err(error, 'defining time units' )

!--- define longitude
error = nf90_def_var(ncid, "lon", NF90_FLOAT, (/dim_id_i/), id_lon)
  call netcdf_err(error, 'defining lon' )
error = nf90_put_att(ncid, id_lon, "long_name", "longitude")
  call netcdf_err(error, 'defining lon long name' )
error =  nf90_put_att(ncid, id_lon, "unit", "degrees_east")
  call netcdf_err(error, 'defining lon units' )  

!--- define latitude
error = nf90_def_var(ncid, "lat", NF90_FLOAT, (/dim_id_j/), id_lat)
  call netcdf_err(error, 'defining lat' )
error = nf90_put_att(ncid,  id_lat, "long_name", "latitude")
  call netcdf_err(error, 'defining lat long name' )
error = nf90_put_att(ncid, id_lat, "units", "degrees_north")
   call netcdf_err(error, 'defining lat units' )

!======================================
! percent of snow

error = nf90_def_var(ncid, 'snow_cover_fraction', NF90_FLOAT, (/dim_id_i,dim_id_j, dim_time/), varid)
 call netcdf_err(error, 'defining snow_cover_fraction' )

error = nf90_put_att(ncid, varid, "long_name", "snow cover fraction")
  call netcdf_err(error, 'defining snow_cover_fraction long name' )

error = nf90_put_att(ncid, varid, "units", "%")
  call netcdf_err(error, 'defining snow_cover_fraction units' )

error =  nf90_put_att(ncid, varid, "missing_value", -9999.0)
 call netcdf_err(error, 'adding snow_cover_fraction filling' )

error = nf90_enddef(ncid)
 call netcdf_err(error, 'defining output_filename' )

! --- put time, lat, lon data
error = nf90_put_var(ncid, id_time, sec_since)
    call netcdf_err(error, 'writing time record')

error =  nf90_put_var(ncid, id_lon, longitude, start = (/1/), count = (/destination_lons/))
    call netcdf_err(error, 'writing lon record')

error =  nf90_put_var(ncid, id_lat, latitude, start = (/1/), count = (/destination_lats/))

  call netcdf_err(error, 'writing lat record')

error =  nf90_put_var(ncid, varid , pt_snow_cover, start = (/1,1,1/), count = (/destination_lons,destination_lats,1/))
 call netcdf_err(error, 'writing imssnow_regrid variable')

error = nf90_close(ncid)

end program

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
 
