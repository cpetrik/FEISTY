#/bin/bash

IFILE=/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19610101-20101231.ocean_cobalt_tracers_month_z_FishMIP_CP.nc
OFILE=/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19610101-20101231.ocean_cobalt_tracers_month_z_FishMIP_CP_remapped.nc
STATIC_FILE=full.ocean_static_FishMIP.nc

cdo -O -merge $IFILE $STATIC_FILE $IFILE.tmp

ncatted -O -h \
	-a coordinates,nlgz,c,c,"geolat geolon" \
	-a coordinates,nmdz,c,c,"geolat geolon" \
	$IFILE.tmp

cdo -selname,nlgz,nmdz -sellonlatbox,-180,180,-90,90 -remapbil,grid.des \
	$IFILE.tmp \
    	$OFILE

