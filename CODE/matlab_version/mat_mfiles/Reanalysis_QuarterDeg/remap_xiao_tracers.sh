#/bin/bash
##SBATCH --time=300

iyear=19610101-20101231
STATIC_FILE=/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/full.ocean_static_FishMIP.nc
IFILE=/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/$iyear.ocean_cobalt_tracers_month_z_FishMIP_CP.nc
OFILE=/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/$iyear.ocean_cobalt_tracers_month_z_FishMIP_CP_remapped.nc
TFILE=/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/$iyear.ocean_cobalt_tracers_month_z_FishMIP_CP_temp.nc
cdo -O -merge $IFILE $STATIC_FILE $TFILE
ncatted -O -a coordinates,nlgz,c,c,"geolat geolon" $TFILE
ncatted -O -a coordinates,nmdz,c,c,"geolat geolon" $TFILE
ncks -h -A -v geolat,geolon $STATIC_FILE $TFILE
cdo remapbil,r1440x720 -selname,nlgz,nmdz $TFILE $OFILE
rm -f $TFILE

