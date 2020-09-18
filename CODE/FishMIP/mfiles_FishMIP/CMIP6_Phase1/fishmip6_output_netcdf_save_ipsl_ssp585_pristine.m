% Calc Fish-MIP outputs saved as NetCDF
% Future time period
% SSP 585

clear all
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';
fpath=['/Volumes/FEISTY/NC/FishMIP/IPSL_CMIP6/' cfile '/'];
harv = 'pristine';

%% Fish-MIP OUTPUTS =================================================

load([fpath 'SSP585_empHP_fishMIP_outputs_monthly_' cfile '.mat'])

% PREFERRED (all units = gWW/m2)
%total pelagic biomass tpb = 360x180xMOs
% allPel
%total demersal biomass tdb = 360x180xMOs
% allD
%total consumber biomass tcb = 360x180xMOs
% allC 

% SECONDARY
%total pelagic (Linf <30cm) biomass bp30cm = 360x180xMOs
% SPel 
%total pelagic (>=30 cm and <90cm) biomass bp30to90cm = 360x180xMOs
% MLPel 
%total pelagic (>=90cm) biomass bp90cm = 360x180xMOs
% MLPel 
%total demersal (Linf <30cm) biomass bd30cm = 360x180xMOs
% SDem 
%total demersal (>=30 cm and <90cm) biomass bd30to90cm = 360x180xMOs
% MLDem 
%total demersal (>=90cm) biomass bd90cm = 360x180xMOs
% MLDem

t=(1:length(time))/12;
t=t+2015;

%% Grid
load('/Volumes/FEISTY/Fish-MIP/CMIP6/IPSL/Data_grid_ipsl.mat','GRD');
load('/Volumes/FEISTY/Fish-MIP/CMIP6/IPSL/gridspec_ipsl_cmip6.mat');

[ni,nj] = size(LAT);

%% Reshape to lat,lon,yr
[nid,nt] = size(allC);

tpb = 1.0e20*ones(ni,nj,nt);
tdb = tpb;
tcb = tpb;
bp30cm = tpb;
% bp30to90cm = tpb;
bp90cm = tpb;
% bd30cm = tpb;
% bd30to90cm = tpb;
bd90cm = tpb;

for y=1:nt
    gtpb = 1.0e20*ones(ni,nj);
    ttpb = allPel(:,y);
    gtpb(GRD.ID) = ttpb;
    tpb(:,:,y) = gtpb;
    
    gtdb = 1.0e20*ones(ni,nj);
    ttdb = allD(:,y);
    gtdb(GRD.ID) = ttdb;
    tdb(:,:,y) = gtdb;
    
    gtcb = 1.0e20*ones(ni,nj);
    ttcb = allC(:,y);
    gtcb(GRD.ID) = ttcb;
    tcb(:,:,y) = gtcb;
    
    gp30cm = 1.0e20*ones(ni,nj);
    tp30cm = SPel(:,y);
    gp30cm(GRD.ID) = tp30cm;
    bp30cm(:,:,y) = gp30cm;
    
    gp90cm = 1.0e20*ones(ni,nj);
    tp90cm = LPel(:,y);
    gp90cm(GRD.ID) = tp90cm;
    bp90cm(:,:,y) = gp90cm;
    
%     gd30cm = 1.0e20*ones(ni,nj);
%     td30cm = SDem(:,y);
%     gd30cm(GRD.ID) = td30cm;
%     bd30cm(:,:,y) = gd30cm;
    
    gd90cm = 1.0e20*ones(ni,nj);
    td90cm = allD(:,y);
    gd90cm(GRD.ID) = td90cm;
    bd90cm(:,:,y) = gd90cm;
end

%% Quick look
pb = tpb(:,:,150);
db = tdb(:,:,150);
cb = tcb(:,:,150);
bp30 = bp30cm(:,:,150);
bp90 = bp90cm(:,:,150);
% bd30 = bd30cm(:,:,150);
bd90 = bd90cm(:,:,150);

figure(1)
pcolor(log10(pb'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('allPel')

figure(2)
pcolor(log10(db'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('allDem')

figure(3)
pcolor(log10(cb'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('all Consumers')

figure(4)
pcolor(log10(bp30'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('S Pel')

figure(5)
pcolor(log10(bp90'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('L Pel')

% figure(6)
% pcolor(log10(bd30'))
% shading flat
% colormap('jet')
% colorbar
% caxis([-2 2])

figure(7)
pcolor(log10(bd90'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('L Dem')


%% Output data naming conventions
%<model>_<climate-forcing>_<bias-adjustment>_<climate-scenario>_
%<soc-scenario>_<sens-scenario>_<variable>_<global>_<timestep>_<start-year>_
%<end-year>.nc

%e.g.
%apecosm_ipsl-esm4_nobc_picontrol_histsoc_default_tcb_global_monthly_2001_2010.nc

%% Setup netcdf path to store to
fname1 = 'feisty_ipsl-cm6a_nobc_ssp585_nat_default_';
fname2 = '_global_monthly_2015-2100.nc';

file_tpb = [fpath fname1 'tpb' fname2];
file_tdb = [fpath fname1 'tdb' fname2];
file_tcb = [fpath fname1 'tcb' fname2];
file_bp30 = [fpath fname1 'bp30cm' fname2];
% file_bp3090 = [fpath fname1 'bp30to90cm' fname2];
file_bp90 = [fpath fname1 'bp90cm' fname2];
% file_bd30 = [fpath fname1 'bd30cm' fname2];
% file_bd3090 = [fpath fname1 'bd30to90cm' fname2];
file_bd90 = [fpath fname1 'bd90cm' fname2];

[ni,nj,nt] = size(tpb);

%% tpb
ncidSB = netcdf.create(file_tpb,'NETCDF4');

lon_dim = netcdf.defDim(ncidSB,'longitude',ni);
lat_dim = netcdf.defDim(ncidSB,'latitude',nj);
time_dim = netcdf.defDim(ncidSB,'time',nt);

vidlat = netcdf.defVar(ncidSB,'LAT','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'units','degrees');

vidlon = netcdf.defVar(ncidSB,'LON','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidSB,vidlon,'units','degrees' );

vidtSB = netcdf.defVar(ncidSB,'time','double',time_dim);
netcdf.putAtt(ncidSB,vidtSB,'long_name','time');
netcdf.putAtt(ncidSB,vidtSB,'units','years' );

vidbioSB = netcdf.defVar(ncidSB,'tpb','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSB,vidbioSB,'long_name','total pelagic biomass density');
netcdf.putAtt(ncidSB,vidbioSB,'units','grams wet weight m-2' );
netcdf.defVarFill(ncidSB,vidbioSB,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidSB);

netcdf.putVar(ncidSB,vidlat,LAT);
netcdf.putVar(ncidSB,vidlon,LON);
netcdf.putVar(ncidSB,vidbioSB,tpb);
netcdf.putVar(ncidSB,vidtSB,t);

netcdf.close(ncidSB);

%% tdb
ncidSD = netcdf.create(file_tdb,'NETCDF4');

lon_dim = netcdf.defDim(ncidSD,'longitude',ni);
lat_dim = netcdf.defDim(ncidSD,'latitude',nj);
time_dim = netcdf.defDim(ncidSD,'time',nt);

vidlat = netcdf.defVar(ncidSD,'LAT','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSD,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSD,vidlat,'units','degrees');

vidlon = netcdf.defVar(ncidSD,'LON','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSD,vidlon,'long_name','longitude');
netcdf.putAtt(ncidSD,vidlon,'units','degrees' );

vidtSD = netcdf.defVar(ncidSD,'time','double',time_dim);
netcdf.putAtt(ncidSD,vidtSD,'long_name','time');
netcdf.putAtt(ncidSD,vidtSD,'units','years' );

vidbioSD = netcdf.defVar(ncidSD,'tdb','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSD,vidbioSD,'long_name','total demersal biomass density');
netcdf.putAtt(ncidSD,vidbioSD,'units','grams wet weight m-2' );
netcdf.defVarFill(ncidSD,vidbioSD,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSD,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSD,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidSD,varid,'includes benthos','no');

netcdf.endDef(ncidSD);

netcdf.putVar(ncidSD,vidlat,LAT);
netcdf.putVar(ncidSD,vidlon,LON);
netcdf.putVar(ncidSD,vidbioSD,tdb);
netcdf.putVar(ncidSD,vidtSD,t);

netcdf.close(ncidSD);

%% tcb
ncidCB = netcdf.create(file_tcb,'NETCDF4');

lon_dim = netcdf.defDim(ncidCB,'longitude',ni);
lat_dim = netcdf.defDim(ncidCB,'latitude',nj);
time_dim = netcdf.defDim(ncidCB,'time',nt);

vidlat = netcdf.defVar(ncidCB,'LAT','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidCB,vidlat,'units','degrees');

vidlon = netcdf.defVar(ncidCB,'LON','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidCB,vidlon,'units','degrees' );

vidtCB = netcdf.defVar(ncidCB,'time','double',time_dim);
netcdf.putAtt(ncidCB,vidtCB,'long_name','time');
netcdf.putAtt(ncidCB,vidtCB,'units','years' );

vidbioCB = netcdf.defVar(ncidCB,'tcb','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidCB,vidbioCB,'long_name','total consumer biomass density');
netcdf.putAtt(ncidCB,vidbioCB,'units','grams wet weight m-2' );
netcdf.defVarFill(ncidCB,vidbioCB,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidCB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidCB,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidCB,varid,'includes benthos','yes');

netcdf.endDef(ncidCB);

netcdf.putVar(ncidCB,vidlat,LAT);
netcdf.putVar(ncidCB,vidlon,LON);
netcdf.putVar(ncidCB,vidbioCB,tcb);
netcdf.putVar(ncidCB,vidtCB,t);

netcdf.close(ncidCB);

%% bp30cm
ncid30 = netcdf.create(file_bp30,'NETCDF4');

lon_dim = netcdf.defDim(ncid30,'longitude',ni);
lat_dim = netcdf.defDim(ncid30,'latitude',nj);
time_dim = netcdf.defDim(ncid30,'time',nt);

vidlat = netcdf.defVar(ncid30,'LAT','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid30,vidlat,'long_name','latitude');
netcdf.putAtt(ncid30,vidlat,'units','degrees');

vidlon = netcdf.defVar(ncid30,'LON','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid30,vidlon,'long_name','longitude');
netcdf.putAtt(ncid30,vidlon,'units','degrees' );

vidt30 = netcdf.defVar(ncid30,'time','double',time_dim);
netcdf.putAtt(ncid30,vidt30,'long_name','time');
netcdf.putAtt(ncid30,vidt30,'units','years' );

vidbio30 = netcdf.defVar(ncid30,'bp30cm','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncid30,vidbio30,'long_name','biomass density of pelagic < 30cm');
netcdf.putAtt(ncid30,vidbio30,'units','grams wet weight  m-2' );
netcdf.defVarFill(ncid30,vidbio30,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid30,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid30,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncid30);

netcdf.putVar(ncid30,vidlat,LAT);
netcdf.putVar(ncid30,vidlon,LON);
netcdf.putVar(ncid30,vidbio30,bp30cm);
netcdf.putVar(ncid30,vidt30,t);

netcdf.close(ncid30);

%% bp90cm
ncid90 = netcdf.create(file_bp90,'NETCDF4');

lon_dim = netcdf.defDim(ncid90,'longitude',ni);
lat_dim = netcdf.defDim(ncid90,'latitude',nj);
time_dim = netcdf.defDim(ncid90,'time',nt);

vidlat = netcdf.defVar(ncid90,'LAT','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid90,vidlat,'long_name','latitude');
netcdf.putAtt(ncid90,vidlat,'units','degrees');

vidlon = netcdf.defVar(ncid90,'LON','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid90,vidlon,'long_name','longitude');
netcdf.putAtt(ncid90,vidlon,'units','degrees' );

vidt90 = netcdf.defVar(ncid90,'time','double',time_dim);
netcdf.putAtt(ncid90,vidt90,'long_name','time');
netcdf.putAtt(ncid90,vidt90,'units','years' );

vidbio90 = netcdf.defVar(ncid90,'bp90cm','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncid90,vidbio90,'long_name','biomass density of pelagic >=90cm');
netcdf.putAtt(ncid90,vidbio90,'units','grams wet weight  m-2' );
netcdf.defVarFill(ncid90,vidbio90,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid90,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid90,varid,'wet weight:C ratio','9:1');
% netcdf.putAtt(ncid90,varid,'feisty fish size','29.24 to 232.08cm');

netcdf.endDef(ncid90);

netcdf.putVar(ncid90,vidlat,LAT);
netcdf.putVar(ncid90,vidlon,LON);
netcdf.putVar(ncid90,vidbio90,bp90cm);
netcdf.putVar(ncid90,vidt90,t);

netcdf.close(ncid90);

%% bd90cm
ncid90 = netcdf.create(file_bd90,'NETCDF4');

lon_dim = netcdf.defDim(ncid90,'longitude',ni);
lat_dim = netcdf.defDim(ncid90,'latitude',nj);
time_dim = netcdf.defDim(ncid90,'time',nt);

vidlat = netcdf.defVar(ncid90,'LAT','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid90,vidlat,'long_name','latitude');
netcdf.putAtt(ncid90,vidlat,'units','degrees');

vidlon = netcdf.defVar(ncid90,'LON','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid90,vidlon,'long_name','longitude');
netcdf.putAtt(ncid90,vidlon,'units','degrees' );

vidt90 = netcdf.defVar(ncid90,'time','double',time_dim);
netcdf.putAtt(ncid90,vidt90,'long_name','time');
netcdf.putAtt(ncid90,vidt90,'units','years' );

vidbio90 = netcdf.defVar(ncid90,'bd90cm','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncid90,vidbio90,'long_name','biomass density of demersal >=90cm');
netcdf.putAtt(ncid90,vidbio90,'units','grams wet weight  m-2' );
netcdf.defVarFill(ncid90,vidbio90,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid90,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid90,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncid90,varid,'includes benthos','no');
% netcdf.putAtt(ncid90,varid,'feisty fish size','29.24 to 232.08cm');

netcdf.endDef(ncid90);

netcdf.putVar(ncid90,vidlat,LAT);
netcdf.putVar(ncid90,vidlon,LON);
netcdf.putVar(ncid90,vidbio90,bd90cm);
netcdf.putVar(ncid90,vidt90,t);

netcdf.close(ncid90);







