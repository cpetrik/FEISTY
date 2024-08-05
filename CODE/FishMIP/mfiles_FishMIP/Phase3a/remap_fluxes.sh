#/bin/bash

IFILE=/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CP.nc
OFILE=/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.nc
STATIC_FILE=full.ocean_static_FishMIP.nc

cdo -O -merge $IFILE $STATIC_FILE $IFILE.tmp

ncatted -O -h \
	-a coordinates,jhploss_nlgz_100,c,c,"geolat geolon" \
	-a coordinates,jhploss_nmdz_100,c,c,"geolat geolon" \
	$IFILE.tmp

cdo -selname,jhploss_nlgz_100,jhploss_nmdz_100 -sellonlatbox,-180,180,-90,90 -remapbil,grid.des \
	$IFILE.tmp \
    	$OFILE

