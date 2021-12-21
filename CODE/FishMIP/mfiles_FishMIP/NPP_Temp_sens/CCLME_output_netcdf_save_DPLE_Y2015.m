% Calc Fish-MIP outputs saved as NetCDF
% Historical time period

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/MIP/FEISTY_output/FishMIP/CESM1-BEC/' cfile '/'];
harv = 'pristine';

%% Fish-MIP OUTPUTS =================================================

load([fpath 'CESM_Hist_2000_2006_empHP_fishMIP_outputs_monthly_' cfile '.mat'])

t=time;
mo=t/12;
mo=mo+1850;
yr20=find(mo>2000);
time = mo(yr20);

%% Grid
cpath = '/Volumes/MIP/Fish-MIP/CMIP5/CESM/';
load([cpath 'gridspec_cesm.mat']);
load([cpath 'Data_grid_cesm.mat']);

[ni,nj] = size(LAT);

%% Reshape to lat,lon,yr
[nid,nt] = size(allB);

%%
tpb = 1.000000020040877e+20*ones(ni,nj,nt);
tdb = tpb;
tfb = tpb;
tbb = tpb;

for y=1:nt
    gtpb = 1.000000020040877e+20*ones(ni,nj);
    ttpb = allP(:,y);
    gtpb(GRD.ID) = ttpb;
    tpb(:,:,y) = gtpb;

    gtdb = 1.000000020040877e+20*ones(ni,nj);
    ttdb = allD(:,y);
    gtdb(GRD.ID) = ttdb;
    tdb(:,:,y) = gtdb;

    gtfb = 1.000000020040877e+20*ones(ni,nj);
    ttfb = allF(:,y);
    gtfb(GRD.ID) = ttfb;
    tfb(:,:,y) = gtfb;

    gtbb = 1.000000020040877e+20*ones(ni,nj);
    ttbb = allB(:,y);
    gtbb(GRD.ID) = ttbb;
    tbb(:,:,y) = gtbb;
end

%% Quick look
pb = tpb(:,:,50);
db = tdb(:,:,50);
cb = tfb(:,:,50);
bp30 = tbb(:,:,50);

figure(1)
pcolor(log10(pb'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('allP')

figure(2)
pcolor(log10(db'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('allD')

figure(3)
pcolor(log10(cb'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('all F')

figure(4)
pcolor(log10(bp30'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('All B')

figure(5)
pcolor(LAT)
shading flat
title('Lat')

figure(6)
pcolor(LON)
shading flat
title('Lon')

%% Output data naming conventions
%<model>_<climate-forcing>_<bias-adjustment>_<climate-scenario>_
%<soc-scenario>_<sens-scenario>_<variable>_<global>_<timestep>_<start-year>_
%<end-year>.nc

%e.g.
%apecosm_ipsl-esm4_nobasd_picontrol_histsoc_default_tcb_global_monthly_2001_2010.nc

close all

%% Setup netcdf path to store to
fname1 = 'feisty_cesm1-bgc_nobc_historical_nosoc_co2_nofishing_';
fname2 = '_global_monthly_2000_2006.nc';

file_tpb = [fpath fname1 'tpb' fname2];
file_tdb = [fpath fname1 'tdb' fname2];
file_tfb = [fpath fname1 'tfb' fname2];
file_tbb = [fpath fname1 'tbb' fname2];

[ni,nj,nt] = size(tpb);

%% 
%Lat & Lon should be vectors
LAT = LAT(1,:);
LON = LON(:,1);

%Use Netcdf4 classic
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

%% tpb
ncidSB = netcdf.create(file_tpb,cmode);

time_dim = netcdf.defDim(ncidSB,'time',nt);
lon_dim = netcdf.defDim(ncidSB,'lon',ni);
lat_dim = netcdf.defDim(ncidSB,'lat',nj);

vidtSB = netcdf.defVar(ncidSB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidSB,vidtSB,'long_name','time');
netcdf.putAtt(ncidSB,vidtSB,'standard_name','time');
netcdf.putAtt(ncidSB,vidtSB,'units','year' );
netcdf.putAtt(ncidSB,vidtSB,'calendar','365_day');
netcdf.putAtt(ncidSB,vidtSB,'axis','T');

vidlon = netcdf.defVar(ncidSB,'lon','NC_DOUBLE',lon_dim);
netcdf.putAtt(ncidSB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidSB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidSB,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidSB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidSB,'lat','NC_DOUBLE',lat_dim);
netcdf.putAtt(ncidSB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidSB,vidlat,'axis','Y');

vidbioSB = netcdf.defVar(ncidSB,'tpb','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidSB,vidbioSB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidSB,vidbioSB,'long_name','Biomass of Large Pelagic Fish');
netcdf.putAtt(ncidSB,vidbioSB,'units','g m-2' );
netcdf.defVarFill(ncidSB,vidbioSB,false,1.000000020040877e+20);
netcdf.putAtt(ncidSB,vidbioSB,'missing value',1.000000020040877e+20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSB,varid,'_FillValue',1.000000020040877e+20);
netcdf.putAtt(ncidSB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidSB,varid,'institution','UCSD');
netcdf.putAtt(ncidSB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidSB);

tpb = single(tpb);
netcdf.putVar(ncidSB,vidlat,LAT);
netcdf.putVar(ncidSB,vidlon,LON);
netcdf.putVar(ncidSB,vidbioSB,tpb);
netcdf.putVar(ncidSB,vidtSB,time);

netcdf.close(ncidSB);
%%
ncdisp(file_tpb)

%% tdb
ncidSD = netcdf.create(file_tdb,cmode);

time_dim = netcdf.defDim(ncidSD,'time',nt);
lon_dim = netcdf.defDim(ncidSD,'lon',ni);
lat_dim = netcdf.defDim(ncidSD,'lat',nj);

vidtSD = netcdf.defVar(ncidSD,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidSD,vidtSD,'long_name','time');
netcdf.putAtt(ncidSD,vidtSD,'standard_name','time');
netcdf.putAtt(ncidSD,vidtSD,'calendar','365_day');
netcdf.putAtt(ncidSD,vidtSD,'axis','T');
netcdf.putAtt(ncidSD,vidtSD,'units','year' );

vidlon = netcdf.defVar(ncidSD,'lon','NC_DOUBLE',lon_dim);
netcdf.putAtt(ncidSD,vidlon,'long_name','longitude');
netcdf.putAtt(ncidSD,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidSD,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidSD,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidSD,'lat','NC_DOUBLE',lat_dim);
netcdf.putAtt(ncidSD,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSD,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidSD,vidlat,'units','degrees_north');
netcdf.putAtt(ncidSD,vidlat,'axis','Y');

vidbioSD = netcdf.defVar(ncidSD,'tdb','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidSD,vidbioSD,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidSD,vidbioSD,'long_name','Biomass of Demersal Fish');
netcdf.putAtt(ncidSD,vidbioSD,'units','g m-2' );
netcdf.defVarFill(ncidSD,vidbioSD,false,1.000000020040877e+20);
netcdf.putAtt(ncidSD,vidbioSD,'missing value',1.000000020040877e+20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSD,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSD,varid,'_FillValue',1.000000020040877e+20);
netcdf.putAtt(ncidSD,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidSD,varid,'institution','UCSD');
netcdf.putAtt(ncidSD,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidSD,varid,'includes benthos','no');

netcdf.endDef(ncidSD);

tdb = single(tdb);
netcdf.putVar(ncidSD,vidlat,LAT);
netcdf.putVar(ncidSD,vidlon,LON);
netcdf.putVar(ncidSD,vidbioSD,tdb);
netcdf.putVar(ncidSD,vidtSD,time);

netcdf.close(ncidSD);

%% tfb
ncidCB = netcdf.create(file_tfb,cmode);

time_dim = netcdf.defDim(ncidCB,'time',nt);
lon_dim = netcdf.defDim(ncidCB,'lon',ni);
lat_dim = netcdf.defDim(ncidCB,'lat',nj);

vidtCB = netcdf.defVar(ncidCB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidCB,vidtCB,'long_name','time');
netcdf.putAtt(ncidCB,vidtCB,'standard_name','time');
netcdf.putAtt(ncidCB,vidtCB,'calendar','365_day');
netcdf.putAtt(ncidCB,vidtCB,'axis','T');
netcdf.putAtt(ncidCB,vidtCB,'units','year' );

vidlon = netcdf.defVar(ncidCB,'lon','NC_DOUBLE',lon_dim);
netcdf.putAtt(ncidCB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidCB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidCB,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidCB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidCB,'lat','NC_DOUBLE',lat_dim);
netcdf.putAtt(ncidCB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidCB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidCB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidCB,vidlat,'axis','Y');

vidbioCB = netcdf.defVar(ncidCB,'tcb','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidCB,vidbioCB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidCB,vidbioCB,'long_name','Biomass of Forage Fish');
netcdf.putAtt(ncidCB,vidbioCB,'units','g m-2' );
netcdf.defVarFill(ncidCB,vidbioCB,false,1.000000020040877e+20);
netcdf.putAtt(ncidCB,vidbioCB,'missing value',1.000000020040877e+20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidCB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidCB,varid,'_FillValue',1.000000020040877e+20);
netcdf.putAtt(ncidCB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidCB,varid,'institution','UCSD');
netcdf.putAtt(ncidCB,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidCB,varid,'includes benthos','yes');

netcdf.endDef(ncidCB);

tfb = single(tfb);
netcdf.putVar(ncidCB,vidlat,LAT);
netcdf.putVar(ncidCB,vidlon,LON);
netcdf.putVar(ncidCB,vidbioCB,tfb);
netcdf.putVar(ncidCB,vidtCB,time);

netcdf.close(ncidCB);

%% tbb
ncid30 = netcdf.create(file_tbb,cmode);

time_dim = netcdf.defDim(ncid30,'time',nt);
lon_dim = netcdf.defDim(ncid30,'lon',ni);
lat_dim = netcdf.defDim(ncid30,'lat',nj);

vidt30 = netcdf.defVar(ncid30,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncid30,vidt30,'long_name','time');
netcdf.putAtt(ncid30,vidt30,'standard_name','time');
netcdf.putAtt(ncid30,vidt30,'calendar','365_day');
netcdf.putAtt(ncid30,vidt30,'axis','T');
netcdf.putAtt(ncid30,vidt30,'units','year' );

vidlon = netcdf.defVar(ncid30,'lon','NC_DOUBLE',lon_dim);
netcdf.putAtt(ncid30,vidlon,'long_name','longitude');
netcdf.putAtt(ncid30,vidlon,'standard_name','longitude');
netcdf.putAtt(ncid30,vidlon,'units','degrees_east' );
netcdf.putAtt(ncid30,vidlon,'axis','X');

vidlat = netcdf.defVar(ncid30,'lat','NC_DOUBLE',lat_dim);
netcdf.putAtt(ncid30,vidlat,'long_name','latitude');
netcdf.putAtt(ncid30,vidlat,'standard_name','latitude');
netcdf.putAtt(ncid30,vidlat,'units','degrees_north');
netcdf.putAtt(ncid30,vidlat,'axis','Y');

vidbio30 = netcdf.defVar(ncid30,'bp30cm','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncid30,vidbio30,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncid30,vidbio30,'long_name','Biomass of Benthic Invertebrates');
netcdf.putAtt(ncid30,vidbio30,'units','g  m-2' );
netcdf.defVarFill(ncid30,vidbio30,false,1.000000020040877e+20);
netcdf.putAtt(ncid30,vidbio30,'missing value',1.000000020040877e+20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid30,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid30,varid,'_FillValue',1.000000020040877e+20);
netcdf.putAtt(ncid30,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncid30,varid,'institution','UCSD');
netcdf.putAtt(ncid30,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncid30);

tbb = single(tbb);
netcdf.putVar(ncid30,vidlat,LAT);
netcdf.putVar(ncid30,vidlon,LON);
netcdf.putVar(ncid30,vidbio30,tbb);
netcdf.putVar(ncid30,vidt30,time);

netcdf.close(ncid30);


