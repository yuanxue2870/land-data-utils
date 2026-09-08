
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
  
  outname = "/scratch4/NCEPDEV/stmp/Yuan.Xue/snodas/sorc/gaussian_mapping/snodas1km_retro_mapping.nc"
   system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
  outfile = addfile(outname, "c")

  corners_file = addfile("/scratch4/NCEPDEV/stmp/Yuan.Xue/snodas/retro/retro_vegsoil.nc", "r")
  retro_lat = corners_file->lat
  retro_lon = corners_file->lon
   
  snodas_file = addfile("/scratch4/NCEPDEV/land/data/evaluation/SNODAS/orig/SNODAS_unmasked_20231231.nc", "r")
  snodas_lat = snodas_file->lat
  snodas_lon = snodas_file->lon
  snodas_lon = where(snodas_lon.lt.0,snodas_lon+360,snodas_lon)

  pt_i_snodas = new((/8192/),short)
   pt_i_snodas!0 = "lon"
  pt_j_snodas = new((/4096/),short)
  pt_j_snodas!0 = "lat"
 
  pt_i_snodas = toshort(ind_nearest_coord(snodas_lon,retro_lon,0))
  pt_j_snodas = toshort(ind_nearest_coord(snodas_lat,retro_lat,0))
    
  pt_j_snodas = pt_j_snodas + toshort(1)
  pt_i_snodas = pt_i_snodas + toshort(1)
  
  outfile->pt_i_snodas = pt_i_snodas
  outfile->pt_j_snodas = pt_j_snodas
  
end 
