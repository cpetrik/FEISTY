#/bin/bash
##SBATCH --time=300

iyear=19590101-20101231
STATIC_FILE=/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/full.ocean_static_FishMIP.nc
IFILE=/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/$iyear.ocean_cobalt_fluxes_int_FishMIP_CPcopy.nc
OFILE=/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/$iyear.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.nc
TFILE=/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/$iyear.ocean_cobalt_fluxes_int_FishMIP_CP_temp.nc
cdo -O -merge $IFILE $STATIC_FILE $TFILE
ncatted -O -a coordinates,jhploss_nlgz_100,c,c,"geolat geolon" $TFILE
ncatted -O -a coordinates,jhploss_nmdz_100,c,c,"geolat geolon" $TFILE
ncks -h -A -v geolat,geolon $STATIC_FILE $TFILE
cdo remapbil,r1440x720 -selname,jhploss_nlgz_100,jhploss_nmdz_100 $TFILE $OFILE
rm -f $TFILE

