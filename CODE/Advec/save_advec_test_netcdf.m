% SAVE 30 yr advection sun AS NETCDF

dpath = '/Volumes/GFDL/NC/AdvectTests/';

bio = csvread([dpath 'bio_2Dadvect_test_global_vel200_dt1hr_j2_30yrs.csv']);
yrs=[1:length(totb)]/365;

ncidUV = netcdf.create([dpath 'bio_2Dadvect_test_global_vel200_dt1hr_j2_30yrs.nc'],'NC_WRITE');
time_dim = netcdf.defDim(ncidUV,'ntime',10950);
cell_dim = netcdf.defDim(ncidUV,'ncell',48111);

varidT = netcdf.defVar(ncidUV,'time','double',time_dim);
varidB = netcdf.defVar(ncidUV,'biomass','double',[time_dim,cell_dim]);
netcdf.endDef(ncidUV);
netcdf.putVar(ncidUV,varidT,yrs);
netcdf.putVar(ncidUV,varidB,bio);
netcdf.close(ncidUV);