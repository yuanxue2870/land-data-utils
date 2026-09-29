program create_fv3_mapping

  use netcdf
  implicit none

  integer            :: tile_dim 
  character*100      :: tile_path 
  character*100      :: orog_path
  character*100      :: out_path 
  character*20       :: otype ! orography filename stub. For atm only, oro_C${RES}, for atm/ocean oro_C${RES}.mx100
  character*10        :: obs_source
  integer             :: source_i_size ! size of input data 
  integer             :: source_j_size

  logical :: file_exists 

! IMS input info
!  integer             :: length ! binary file recl
  character*100      :: ims_path
  character*100      :: ims_lat_name
  character*100      :: ims_lon_name

  logical            :: include_source_latlon = .false.
  real, parameter    :: perturb_value         = 1.d-4    ! a small adjustment to lat/lon to find [radians]
  integer            :: tile_length
  integer :: fv3_search_order(6) = (/3,1,2,5,6,4/)
  integer :: quick_search_pad = 1

  real   , allocatable, dimension(:,:,:)          :: fv3_lat, fv3_lon 
  real   , allocatable, dimension(:,:,:)          :: fv3_lat_cnt, fv3_lon_cnt, fv3_oro
  integer, allocatable, dimension(:,:)            :: fv3_mask

  integer, allocatable, dimension(:,:) :: lookup_tile, lookup_i, lookup_j
  real   , allocatable, dimension(:,:) :: source_lat, source_lon, source_data
  
  real, dimension(4) :: lat_vertex, lon_vertex
  
  integer :: itile, tile_index, tile_i_index, tile_j_index, source_i_index, source_j_index
  integer :: tile_save, tile_i_save, tile_j_save, pad_i_min, pad_i_max, pad_j_min, pad_j_max
  logical :: found, inside_a_polygon
  real    :: lat2find, lon2find
  integer :: ncid, dimid, varid, ierr   ! netcdf identifiers
  integer :: dim_id_i, dim_id_j           ! netcdf dimension identifiers
  integer :: dim_id_i_fv3, dim_id_j_fv3, dim_id_t_fv3
  integer :: i,j,t, io
  character*250 :: filename
  character*9 :: tilestr
  character*20 :: dimstr
  real, parameter :: deg2rad = 3.1415926535897931/180.0

  namelist/fv3_mapping_nml/ tile_dim, otype, orog_path, tile_path, obs_source, source_i_size, source_j_size, ims_path, ims_lat_name, ims_lon_name, out_path

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Setup inputs and read namelist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! read namelist
 
! defaults
 
 inquire(file='fv3_mapping.nml', exist=file_exists)

 if (.not. file_exists) then
        print *, 'namelistfile does not exist, exiting'
        stop 10
 endif

 open (action='read', file='fv3_mapping.nml', iostat=ierr, newunit=io)
 read (nml=fv3_mapping_nml, iostat=ierr, unit=io)
 close (io)

! write(*,*) tile_dim, otype, orog_path, tile_path, obs_source, source_i_size, source_j_size, ims_path, ims_lat_name, ims_lon_name, out_path

 tile_length= tile_dim*2 + 1

 allocate(fv3_lat(tile_length,tile_length,6))
 allocate(fv3_lon(tile_length,tile_length,6))
 allocate(fv3_lon_cnt(tile_dim,tile_dim,6))
 allocate(fv3_lat_cnt(tile_dim,tile_dim,6))
 allocate(fv3_oro(tile_dim,tile_dim,6))
 allocate(fv3_mask(tile_dim,tile_dim))

 allocate(lookup_tile(source_i_size,source_j_size))  
 allocate(lookup_i(source_i_size,source_j_size))
 allocate(lookup_j(source_i_size,source_j_size))
 allocate(source_lat(source_i_size,source_j_size))
 allocate(source_lon(source_i_size,source_j_size))
 allocate(source_data(source_i_size,source_j_size))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read IMS lat/lon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  if ( obs_source(1:3)=="IMS" ) then 
      write(6,*) 'Reading in IMS coordinate info' 
!      length = source_i_size*source_j_size*8  + 1
      filename = trim(ims_path)//trim(ims_lat_name)
      open(10, file=filename, form='unformatted', access='stream')
      read(10) source_data
      source_lat = -9999.
      where(abs(source_data) <= 90.) source_lat = source_data
      
      filename = trim(ims_path)//trim(ims_lon_name)
      open(10, file=filename, form='unformatted', access='stream')
      read(10) source_data
      source_lon = -9999.
      where(abs(source_data) <= 360.) source_lon = source_data
      where(source_lon < 0.0 .and. source_lon >= -180.0) source_lon = source_lon + 360.0
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read FV3 tile information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(dimstr,*)  tile_dim

  do itile = 1, 6
 
    write(tilestr,'(a5,i1,a3)')  ".tile", itile, ".nc"

    filename = trim(tile_path)//"/C"//trim(adjustl(dimstr))//"_grid"//tilestr
    write(6,*) 'Reading in tile file' , filename

    ierr = nf90_open(filename, NF90_NOWRITE, ncid)
      if (ierr /= nf90_noerr) call handle_err(ierr)

    ierr = nf90_inq_varid(ncid, "x", varid)
    ierr = nf90_get_var(ncid, varid , fv3_lon(:,:,itile))
  
    ierr = nf90_inq_varid(ncid, "y", varid)
    ierr = nf90_get_var(ncid, varid , fv3_lat(:,:,itile))

    ierr = nf90_close(ncid)

    ! get orography

    filename = trim(orog_path)//"/"//trim(otype)//tilestr
    write(6,*) 'Reading in orog file' , filename

    ierr = nf90_open(filename, NF90_NOWRITE, ncid)
      if (ierr /= nf90_noerr) call handle_err(ierr)

    ierr = nf90_inq_varid(ncid, "orog_filt", varid)
    ierr = nf90_get_var(ncid, varid , fv3_oro(:,:,itile))

    ierr = nf90_close(ncid)

  end do
  
! get center of grid cell for output

  do t=1,6
     do i =1,tile_dim 
        do j=1,tile_dim 
                fv3_lon_cnt(i,j,t) = fv3_lon(i*2,j*2,t)
                fv3_lat_cnt(i,j,t) = fv3_lat(i*2,j*2,t)
        enddo 
     enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop through the source points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  lookup_tile = -9999
  lookup_i = -9999
  lookup_j = -9999
  
  tile_save = fv3_search_order(1)
  tile_i_save = 1
  tile_j_save = 1
  
  source_i_loop : do source_i_index = 1, source_i_size
  source_j_loop : do source_j_index = 1, source_j_size
  
    found = .false.
    lat2find = source_lat(source_i_index, source_j_index)
    lon2find = source_lon(source_i_index, source_j_index)
    
    if(lat2find < -90. .or. lat2find > 90.  .or. &
       lon2find <   0. .or. lon2find > 360.) cycle source_j_loop     ! skip if out of projections
   
    ! input is in degrees. 
    lat2find = deg2rad * lat2find
    lon2find = deg2rad * lon2find
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! check around the last found tile/i/j
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    pad_i_min = max(tile_i_save-quick_search_pad,1)
    pad_i_max = min(tile_i_save+quick_search_pad,tile_dim)
    pad_j_min = max(tile_j_save-quick_search_pad,1)
    pad_j_max = min(tile_j_save+quick_search_pad,tile_dim)
    
    tile_index = tile_save
    
    do tile_i_index = pad_i_min, pad_i_max
    do tile_j_index = pad_j_min, pad_j_max
      
      lat_vertex(1) = fv3_lat((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 1,tile_index)  ! LL
      lat_vertex(2) = fv3_lat((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 1,tile_index)  ! LR
      lat_vertex(3) = fv3_lat((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 3,tile_index)  ! UR
      lat_vertex(4) = fv3_lat((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 3,tile_index)  ! UL
      
      lon_vertex(1) = fv3_lon((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 1,tile_index)  ! LL
      lon_vertex(2) = fv3_lon((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 1,tile_index)  ! LR
      lon_vertex(3) = fv3_lon((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 3,tile_index)  ! UR
      lon_vertex(4) = fv3_lon((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 3,tile_index)  ! UL
        
      lat_vertex = lat_vertex * deg2rad
      lon_vertex = lon_vertex * deg2rad
      
      found = inside_a_polygon(lon2find, lat2find, 4, lon_vertex, lat_vertex)
        
      if(found) then
        lookup_tile(source_i_index,source_j_index) = tile_index
        lookup_i   (source_i_index,source_j_index) = tile_i_index
        lookup_j   (source_i_index,source_j_index) = tile_j_index
        tile_save = tile_index
        tile_i_save = tile_i_index
        tile_j_save = tile_j_index
        cycle source_j_loop
      end if
        
    end do
    end do
      
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! not found so do a general check
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do itile = 1, 6

      tile_index = fv3_search_order(itile)
      
      do tile_i_index = 1, tile_dim
      do tile_j_index = 1, tile_dim
      
        lat_vertex(1) = fv3_lat((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 1,tile_index)  ! LL
        lat_vertex(2) = fv3_lat((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 1,tile_index)  ! LR
        lat_vertex(3) = fv3_lat((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 3,tile_index)  ! UR
        lat_vertex(4) = fv3_lat((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 3,tile_index)  ! UL
      
        lon_vertex(1) = fv3_lon((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 1,tile_index)  ! LL
        lon_vertex(2) = fv3_lon((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 1,tile_index)  ! LR
        lon_vertex(3) = fv3_lon((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 3,tile_index)  ! UR
        lon_vertex(4) = fv3_lon((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 3,tile_index)  ! UL
        
        lat_vertex = lat_vertex * deg2rad
        lon_vertex = lon_vertex * deg2rad
      
        found = inside_a_polygon(lon2find, lat2find, 4, lon_vertex, lat_vertex)
        
        if(found) then
          lookup_tile(source_i_index,source_j_index) = tile_index
          lookup_i   (source_i_index,source_j_index) = tile_i_index
          lookup_j   (source_i_index,source_j_index) = tile_j_index
          tile_save = tile_index
          tile_i_save = tile_i_index
          tile_j_save = tile_j_index
          cycle source_j_loop
        end if
        
      end do
      end do
      
    end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! not found so do a general check with a perturbation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    print*, "Did not find, add perturbation"
  
    lat2find = lat2find + perturb_value
    lon2find = lon2find + perturb_value
  
    do itile = 1, 6

      tile_index = fv3_search_order(itile)

      do tile_i_index = 1, tile_dim
      do tile_j_index = 1, tile_dim

        lat_vertex(1) = fv3_lat((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 1,tile_index)  ! LL
        lat_vertex(2) = fv3_lat((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 1,tile_index)  ! LR
        lat_vertex(3) = fv3_lat((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 3,tile_index)  ! UR
        lat_vertex(4) = fv3_lat((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 3,tile_index)  ! UL

        lon_vertex(1) = fv3_lon((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 1,tile_index)  ! LL
        lon_vertex(2) = fv3_lon((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 1,tile_index)  ! LR
        lon_vertex(3) = fv3_lon((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 3,tile_index)  ! UR
        lon_vertex(4) = fv3_lon((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 3,tile_index)  ! UL

        lat_vertex = lat_vertex * deg2rad
        lon_vertex = lon_vertex * deg2rad

        found = inside_a_polygon(lon2find, lat2find, 4, lon_vertex, lat_vertex)
  
        if(found) then
          lookup_tile(source_i_index,source_j_index) = tile_index
          lookup_i   (source_i_index,source_j_index) = tile_i_index
          lookup_j   (source_i_index,source_j_index) = tile_j_index
          tile_save = tile_index
          tile_i_save = tile_i_index
          tile_j_save = tile_j_index
          cycle source_j_loop
        end if

      end do
      end do

    end do
   
    if(.not.found) then
      print*, "Did not find in cube sphere:", source_lat(source_i_index, source_j_index), ",", source_lon(source_i_index, source_j_index)
      stop
    end if

  end do source_j_loop
     if(mod(source_i_index,10) == 0) print *, "finished loop: ",source_i_index, " of ", source_i_size
  end do source_i_loop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create the output filename and netcdf file (overwrite old)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  filename= trim(out_path)//trim(obs_source)//"_to_FV3_mapping."//trim(otype)//".nc"
  write(6,*) 'writing indexes to ', trim(filename)

  ierr = nf90_create(filename, NF90_NETCDF4, ncid)
    if (ierr /= nf90_noerr) call handle_err(ierr)

! Define dimensions in the file.

  ierr = nf90_def_dim(ncid, "idim"   , source_i_size , dim_id_i)
    if (ierr /= nf90_noerr) call handle_err(ierr)
  ierr = nf90_def_dim(ncid, "jdim"   , source_j_size , dim_id_j)
    if (ierr /= nf90_noerr) call handle_err(ierr)

! fv3 lat/lon
  ierr = nf90_def_dim(ncid, "idim_fv3"   , tile_dim , dim_id_i_fv3)
    if (ierr /= nf90_noerr) call handle_err(ierr)
  ierr = nf90_def_dim(ncid, "jdim_fv3"   , tile_dim , dim_id_j_fv3)
    if (ierr /= nf90_noerr) call handle_err(ierr)
  ierr = nf90_def_dim(ncid, "tdim_fv3"   , 6 , dim_id_t_fv3)
    if (ierr /= nf90_noerr) call handle_err(ierr)
  
! Define variables in the file.

  ierr = nf90_def_var(ncid, "tile", NF90_INT, (/dim_id_j, dim_id_i/), varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)

    ierr = nf90_put_att(ncid, varid, "long_name", "fv3 tile location")
      if (ierr /= nf90_noerr) call handle_err(ierr)
    ierr = nf90_put_att(ncid, varid, "missing_value", -9999)
      if (ierr /= nf90_noerr) call handle_err(ierr)

  ierr = nf90_def_var(ncid, "tile_i", NF90_INT, (/dim_id_j, dim_id_i/), varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)

    ierr = nf90_put_att(ncid, varid, "long_name", "fv3 i location in tile")
      if (ierr /= nf90_noerr) call handle_err(ierr)
    ierr = nf90_put_att(ncid, varid, "missing_value", -9999)
      if (ierr /= nf90_noerr) call handle_err(ierr)

  ierr = nf90_def_var(ncid, "tile_j", NF90_INT, (/dim_id_j, dim_id_i/), varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)

    ierr = nf90_put_att(ncid, varid, "long_name", "fv3 j location in tile")
      if (ierr /= nf90_noerr) call handle_err(ierr)
    ierr = nf90_put_att(ncid, varid, "missing_value", -9999)
      if (ierr /= nf90_noerr) call handle_err(ierr)

  ierr = nf90_def_var(ncid, "lon_fv3", NF90_FLOAT, (/dim_id_j_fv3, dim_id_i_fv3,dim_id_t_fv3/), varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)

    ierr = nf90_put_att(ncid, varid, "long_name", "longitude fv3 grid")
      if (ierr /= nf90_noerr) call handle_err(ierr)
    ierr = nf90_put_att(ncid, varid, "missing_value", -9999)
      if (ierr /= nf90_noerr) call handle_err(ierr)

  ierr = nf90_def_var(ncid, "lat_fv3", NF90_FLOAT, (/dim_id_j_fv3, dim_id_i_fv3,dim_id_t_fv3/), varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)

    ierr = nf90_put_att(ncid, varid, "long_name", "latitude fv3 grid")
      if (ierr /= nf90_noerr) call handle_err(ierr)
    ierr = nf90_put_att(ncid, varid, "missing_value", -9999)
      if (ierr /= nf90_noerr) call handle_err(ierr)

  ierr = nf90_def_var(ncid, "oro_fv3", NF90_FLOAT, (/dim_id_j_fv3, dim_id_i_fv3,dim_id_t_fv3/), varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)

    ierr = nf90_put_att(ncid, varid, "long_name", "orography fv3 grid")
      if (ierr /= nf90_noerr) call handle_err(ierr)
    ierr = nf90_put_att(ncid, varid, "missing_value", -9999)
      if (ierr /= nf90_noerr) call handle_err(ierr)

 if(include_source_latlon) then

  ierr = nf90_def_var(ncid, "ims_lat", NF90_FLOAT, (/dim_id_j, dim_id_i/), varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)

    ierr = nf90_put_att(ncid, varid, "long_name", "ims latitude")
      if (ierr /= nf90_noerr) call handle_err(ierr)
    ierr = nf90_put_att(ncid, varid, "missing_value", -9999.)
      if (ierr /= nf90_noerr) call handle_err(ierr)

  ierr = nf90_def_var(ncid, "ims_lon", NF90_FLOAT, (/dim_id_j, dim_id_i/), varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)

    ierr = nf90_put_att(ncid, varid, "long_name", "ims longitude")
      if (ierr /= nf90_noerr) call handle_err(ierr)
    ierr = nf90_put_att(ncid, varid, "missing_value", -9999.)
      if (ierr /= nf90_noerr) call handle_err(ierr)
  
 end if

  ierr = nf90_enddef(ncid)

  ierr = nf90_inq_varid(ncid, "tile", varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)
  ierr = nf90_put_var(ncid, varid , lookup_tile)
    if (ierr /= nf90_noerr) call handle_err(ierr)

  ierr = nf90_inq_varid(ncid, "tile_i", varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)
  ierr = nf90_put_var(ncid, varid , lookup_i)
    if (ierr /= nf90_noerr) call handle_err(ierr)

  ierr = nf90_inq_varid(ncid, "tile_j", varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)
  ierr = nf90_put_var(ncid, varid , lookup_j)
    if (ierr /= nf90_noerr) call handle_err(ierr)

  ierr = nf90_inq_varid(ncid, "lon_fv3", varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)
  ierr = nf90_put_var(ncid, varid , fv3_lon_cnt)
    if (ierr /= nf90_noerr) call handle_err(ierr)

  ierr = nf90_inq_varid(ncid, "lat_fv3", varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)
  ierr = nf90_put_var(ncid, varid , fv3_lat_cnt)
    if (ierr /= nf90_noerr) call handle_err(ierr)

  ierr = nf90_inq_varid(ncid, "oro_fv3", varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)
  ierr = nf90_put_var(ncid, varid , fv3_oro)
    if (ierr /= nf90_noerr) call handle_err(ierr)

 if(include_source_latlon) then

  ierr = nf90_inq_varid(ncid, "ims_lat", varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)
  ierr = nf90_put_var(ncid, varid , source_lat)
    if (ierr /= nf90_noerr) call handle_err(ierr)

  ierr = nf90_inq_varid(ncid, "ims_lon", varid)
    if (ierr /= nf90_noerr) call handle_err(ierr)
  ierr = nf90_put_var(ncid, varid , source_lon)
    if (ierr /= nf90_noerr) call handle_err(ierr)
  
 end if

 ierr = nf90_close(ncid)

end program

  subroutine latlon2xyz(siz,lon, lat, x, y, z)
  implicit none
  integer, intent(in) :: siz
  real, intent(in) :: lon(siz), lat(siz)
  real, intent(out) :: x(siz), y(siz), z(siz)
  
  integer n

  do n = 1, siz
    x(n) = cos(lat(n))*cos(lon(n))
    y(n) = cos(lat(n))*sin(lon(n))
    z(n) = sin(lat(n))
  enddo

  end subroutine latlon2xyz

  FUNCTION spherical_angle(v1, v2, v3)
    implicit none
    real, parameter :: EPSLN30 = 1.e-30
    real, parameter :: PI=3.1415926535897931
    real v1(3), v2(3), v3(3)
    real  spherical_angle
 
    real px, py, pz, qx, qy, qz, ddd;
  
    ! vector product between v1 and v2 
    px = v1(2)*v2(3) - v1(3)*v2(2)
    py = v1(3)*v2(1) - v1(1)*v2(3)
    pz = v1(1)*v2(2) - v1(2)*v2(1)
    ! vector product between v1 and v3 
    qx = v1(2)*v3(3) - v1(3)*v3(2);
    qy = v1(3)*v3(1) - v1(1)*v3(3);
    qz = v1(1)*v3(2) - v1(2)*v3(1);

    ddd = (px*px+py*py+pz*pz)*(qx*qx+qy*qy+qz*qz);
    if ( ddd <= 0.0 ) then
      spherical_angle = 0. 
    else 
      ddd = (px*qx+py*qy+pz*qz) / sqrt(ddd);
      if( abs(ddd-1) < EPSLN30 ) ddd = 1;
      if( abs(ddd+1) < EPSLN30 ) ddd = -1;
      if ( ddd>1. .or. ddd<-1. ) then
    !FIX to correctly handle co-linear points (angle near pi or 0) */
    if (ddd < 0.) then
      spherical_angle = PI
    else
      spherical_angle = 0.
    endif
      else
    spherical_angle = acos( ddd )
      endif
    endif  

    return

  END FUNCTION spherical_angle

  FUNCTION inside_a_polygon(lon1, lat1, npts, lon2, lat2)
    implicit none
    real, parameter :: EPSLN10 = 1.e-10
    real, parameter :: EPSLN8 = 1.e-8
    real, parameter :: PI=3.1415926535897931
    real, parameter :: RANGE_CHECK_CRITERIA=0.05
    real :: anglesum, angle, spherical_angle
    integer i, ip1
    real lon1, lat1
    integer npts
    real lon2(npts), lat2(npts)
    real x2(npts), y2(npts), z2(npts)
    real lon1_1d(1), lat1_1d(1)
    real x1(1), y1(1), z1(1)
    real pnt0(3),pnt1(3),pnt2(3)
    logical inside_a_polygon
    real max_x2,min_x2,max_y2,min_y2,max_z2,min_z2
    !first convert to cartesian grid */
    call latlon2xyz(npts,lon2, lat2, x2, y2, z2);
    lon1_1d(1) = lon1
    lat1_1d(1) = lat1
    call latlon2xyz(1,lon1_1d, lat1_1d, x1, y1, z1);
    inside_a_polygon = .false.
    max_x2 = maxval(x2)
    if( x1(1) > max_x2+RANGE_CHECK_CRITERIA ) return
    min_x2 = minval(x2)
    if( x1(1)+RANGE_CHECK_CRITERIA < min_x2 ) return
    max_y2 = maxval(y2)
    if( y1(1) > max_y2+RANGE_CHECK_CRITERIA ) return
    min_y2 = minval(y2)
    if( y1(1)+RANGE_CHECK_CRITERIA < min_y2 ) return
    max_z2 = maxval(z2)
    if( z1(1) > max_z2+RANGE_CHECK_CRITERIA ) return
    min_z2 = minval(z2)
    if( z1(1)+RANGE_CHECK_CRITERIA < min_z2 ) return

    pnt0(1) = x1(1)
    pnt0(2) = y1(1)
    pnt0(3) = z1(1)
    
    anglesum = 0;
    do i = 1, npts
       if(abs(x1(1)-x2(i)) < EPSLN10 .and. &
          abs(y1(1)-y2(i)) < EPSLN10 .and. &
          abs(z1(1)-z2(i)) < EPSLN10 ) then ! same as the corner point
      inside_a_polygon = .true.
      return
       endif
       ip1 = i+1
       if(ip1>npts) ip1 = 1
       pnt1(1) = x2(i)
       pnt1(2) = y2(i)
       pnt1(3) = z2(i)
       pnt2(1) = x2(ip1)
       pnt2(2) = y2(ip1)
       pnt2(3) = z2(ip1)

       angle = spherical_angle(pnt0, pnt2, pnt1);
!       anglesum = anglesum + spherical_angle(pnt0, pnt2, pnt1);
       anglesum = anglesum + angle
    enddo

    if(abs(anglesum-2*PI) < EPSLN8) then
       inside_a_polygon = .true.
    else
       inside_a_polygon = .false.
    endif

    return
    
  end FUNCTION inside_a_polygon

  subroutine handle_err(ierr)
    use netcdf
    integer, intent ( in) :: ierr
 
    if(ierr /= nf90_noerr) then
      print *, trim(nf90_strerror(ierr))
      stop "Stopped"
    end if
  end subroutine handle_err
