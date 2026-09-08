program data_regrid_vector

  use netcdf
  implicit none

  double precision      :: sec_since

  character*256         :: source_path
  character*256         :: destination_path
  
  character*256         :: source_filename
  character*256         :: destination_filename
  character*256         :: weights_filename

  character*19     :: current_date  ! current date
  character*19     :: since_date = "1970-01-01 00:00:00"
  character*4      :: yyyy, year
  character*2      :: mm, dd
  character*10     :: fv3_grid 

!  integer, parameter           :: days                  = 122
  integer, parameter           :: days                  = 151
  integer, parameter           :: source_lats           = 1801
  integer, parameter           :: source_lons           = 3600
  integer :: destination_locs, weight_locs

  integer :: latloc, lonloc, iwt, ndays
  integer :: i, j, itotal, id, iyear, im, io, ierr
  integer :: offset_ss

  integer, dimension(2)                       :: start, count
  character*10     :: idate(days)
  real   , dimension(source_lons,source_lats,days) :: snowdensity
  real   , dimension(source_lons,source_lats,days) :: swe1
  real   , dimension(source_lons,source_lats) :: swe
  real   , dimension(source_lons,source_lats) :: source_input

  integer,  allocatable :: source_lookup(:), destination_lookup(:)
  real*8 ,  allocatable :: weights(:)
  real   ,  allocatable :: era5SWE(:)
  real   ,  allocatable :: era5SWEWT(:)
  real   ,  allocatable :: era5SnowDensity(:)
  real   ,  allocatable :: era5SnowDensityWT(:)

  real   ,  allocatable :: xcb(:)    ! destination lon
  real   ,  allocatable :: ycb(:)    ! destination lat

  integer :: ncid, dimid, varid, status   ! netcdf identifiers
  integer :: dim_id_i, dim_id_t           ! netcdf dimension identifiers

  logical :: file_exists

  real*4  fillVal;
  fillVal = -9999.0

  offset_ss=0

  namelist/regrid_era5_nml/ destination_locs, weight_locs, fv3_grid, source_path, weights_filename, destination_path

! read namelist

  inquire(file='regrid_era5land.nml', exist=file_exists)

  if (.not. file_exists) then
      print *, 'namelistfile does not exist, exiting'
      stop 10
  endif

  open (action='read', file='regrid_era5land.nml', iostat=ierr, newunit=io)
  read (nml=regrid_era5_nml, iostat=ierr, unit=io)
  close (io)

! read input file for dates
  open(18, file='date_input.txt', status='old')
  do id=1,days
    read(18,'(A8)') idate(id)
  end do
  close(18)

  allocate(source_lookup(weight_locs))
  allocate(destination_lookup(weight_locs))
  allocate(weights(weight_locs))
  allocate(era5SWE(destination_locs))
  allocate(era5SWEWT(destination_locs))
  allocate(era5SnowDensity(destination_locs))
  allocate(era5SnowDensityWT(destination_locs))
  allocate(xcb(destination_locs))
  allocate(ycb(destination_locs))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read weights file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  write(*,*) trim(weights_filename)
  status = nf90_open(trim(weights_filename), NF90_NOWRITE, ncid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_inq_varid(ncid, "col", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, varid , source_lookup)
    if (status /= nf90_noerr) call handle_err(status)
  
  status = nf90_inq_varid(ncid, "row", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, varid , destination_lookup)
    if (status /= nf90_noerr) call handle_err(status)
  
  status = nf90_inq_varid(ncid, "S", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, varid , weights)
    if (status /= nf90_noerr) call handle_err(status)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read lattitude and longitude data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   status = nf90_inq_varid(ncid, "yc_b", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, varid , ycb)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "xc_b", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, varid , xcb)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_close(ncid)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read source data in NetCDF format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  source_filename=trim(source_path)//'ERA5land_00Z_1Jan_31May26.nc'
!  write(*,*) trim(source_filename)
  status = nf90_open(trim(source_filename), NF90_NOWRITE, ncid)
    if (status /= nf90_noerr) call handle_err(status)

! read snow water equivalent
  status = nf90_inq_varid(ncid, "sd", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, varid , swe1)
    if (status /= nf90_noerr) call handle_err(status)
   
! read snow density

  status = nf90_inq_varid(ncid, "rsn", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, varid , snowdensity)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_close(ncid)

  do id = 1, days 
 
  destination_filename=trim(destination_path)//trim(fv3_grid)//'_era5land_snow_'//trim(idate(id))//'00.nc'
  
  yyyy=idate(id)(1:4)
  mm=idate(id)(5:6)
  dd=idate(id)(7:8)

  write(*,*) yyyy, mm, dd
  current_date=yyyy//'-'//mm//'-'//dd//' 00:00:00'
  call calc_sec_since(since_date, current_date, offset_ss, sec_since)

!!!!!!!!!! Initialize swe and source_input !!!!!!!!!!!!!!!!!!!!!!!!!!
  do j=1,1801
   do i=1,3600
    swe(i,j)=-9999.0
    source_input(i,j)=-9999.0
   end do
  end do

  do j=1,1801
   do i=1,3600
     swe(i,j)=swe1(i,j,id)
     source_input(i,j)=snowdensity(i,j,id)
   end do
  end do
  !write(*,*) source_input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Regrid the data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  era5SnowDensity = 0.0
  era5SnowdensityWT = 0.0
  era5SWE=0.0
  era5SWEWT=0.0
  
  do iwt = 1, weight_locs

    latloc = source_lookup(iwt)/source_lons + 1
    lonloc = source_lookup(iwt) - (latloc-1)*source_lons
       
    if(source_input(lonloc,latloc) > 0.0 ) then
     era5SnowDensity(destination_lookup(iwt)) =  era5SnowDensity(destination_lookup(iwt)) + weights(iwt) * source_input(lonloc,latloc)
     era5SnowdensityWT(destination_lookup(iwt)) = era5SnowdensityWT(destination_lookup(iwt)) + weights(iwt)
    end if

    if(swe(lonloc,latloc) >= 0.0 ) then
      era5SWE(destination_lookup(iwt)) = era5SWE(destination_lookup(iwt)) + weights(iwt) * swe(lonloc,latloc)
      era5SWEWT(destination_lookup(iwt)) = era5SWEWT(destination_lookup(iwt)) + weights(iwt)
    end if

  end do
  
  do iwt = 1, destination_locs

    if(era5SWEWT(iwt) > 0.0 ) then
     era5SWE(iwt) = era5SWE(iwt)/era5SWEWT(iwt)
    else
     era5SWE(iwt) = -9999.0
    end if

   if(era5SnowdensityWT(iwt) > 0.0 ) then
      era5SnowDensity(iwt) =  era5SnowDensity(iwt)/era5SnowdensityWT(iwt)
    else
      era5SnowDensity(iwt) = -9999.0
    end if

  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write the destination file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  status = nf90_create(destination_filename, NF90_CLOBBER, ncid)
    if (status /= nf90_noerr) call handle_err(status)

! Define dimensions in the file.
 
  status = nf90_def_dim(ncid, "location", destination_locs, dim_id_i)
    if (status /= nf90_noerr) call handle_err(status)  
  status = nf90_def_dim(ncid, "time", NF90_UNLIMITED , dim_id_t)
    if (status /= nf90_noerr) call handle_err(status)
  
! Define variables in the file.

  status = nf90_def_var(ncid, "time", NF90_DOUBLE, dim_id_t, varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "long_name", "time")
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "units", "seconds since "//since_date)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "lat", NF90_FLOAT, (/dim_id_i/), varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "long_name", "latitude")
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "units", "degrees_north")
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "lon", NF90_FLOAT, (/dim_id_i/), varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "long_name", "longitude")
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "unit", "degrees_east")
    if (status /= nf90_noerr) call handle_err(status) 

  status = nf90_def_var(ncid, "SWE", NF90_FLOAT, (/dim_id_i, dim_id_t/), varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "long_name", "era5-land snow water equivalent")
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "units", "m")
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "_FillValue", fillVal)
    if (status /= nf90_noerr) call handle_err(status)

   status = nf90_def_var(ncid, "SnowDensity", NF90_FLOAT, (/dim_id_i, dim_id_t/), varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "long_name", "era5-land snow density")
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "units", "kg/m3")
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "_FillValue", fillVal)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_enddef(ncid)
    if (status /= nf90_noerr) call handle_err(status)

! Write variables in the file.
  
  status = nf90_inq_varid(ncid, "lat", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , ycb, start = (/1/), count = (/destination_locs/))
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "lon", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , xcb, start = (/1/), count = (/destination_locs/))
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "time", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , sec_since)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "SWE", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , era5SWE, start = (/1,1/), count = (/destination_locs,1/))
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "_FillValue", fillVal)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "SnowDensity", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , era5SnowDensity, start = (/1,1/), count = (/destination_locs,1/))
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "_FillValue", fillVal)
    if (status /= nf90_noerr) call handle_err(status)
 
  status = nf90_close(ncid)
    
  enddo                ! end of ndays loop

end program

  subroutine handle_err(status)
    use netcdf
    integer, intent ( in) :: status
 
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine handle_err

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
