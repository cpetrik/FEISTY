% Calc Fish-MIP outputs saved as NetCDF
% PreIndust control

clear all
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';
fpath=['/Volumes/MIP/NC/FishMIP/IPSL_CMIP6/' cfile '/'];
harv = 'pristine';

%% Fish-MIP OUTPUTS =================================================

load([fpath 'PreIndust_empHP_fishMIP_outputs_monthly_' cfile '.mat'])

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
t=t+1950;
clear time

%% Grid
load('/Volumes/MIP/Fish-MIP/CMIP6/IPSL/Data_grid_ipsl.mat','GRD');
load('/Volumes/MIP/Fish-MIP/CMIP6/IPSL/gridspec_ipsl_cmip6.mat');

[ni,nj] = size(LAT);

load('/Volumes/MIP/Fish-MIP/CMIP6/IPSL/preindust/ipsl_pi_temp_btm_monthly_1950_2100.mat',...
    'yr','time');
t_all = time;

runs = find(yr>1950 & yr<=2100);
year = yr(runs);
time = time(runs);
runs1 = find(year>1950 & year<=2015);
runs2 = find(year>2015 & year<=2100);
htime = time(runs1);
ftime = time(runs2);

length(runs1) +length(runs2)

%% Just hsitoric years
allPel = allPel(:,runs1);
allD = allD(:,runs1);
allC = allC(:,runs1);
SPel = SPel(:,runs1);
LPel = LPel(:,runs1);

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
%apecosm_ipsl-esm4_nobasd_picontrol_histsoc_default_tcb_global_monthly_2001_2010.nc

%% Setup netcdf path to store to
fname1 = 'feisty_ipsl-cm6a-lr_nobasd_picontrol_nat_default_';
fname2 = '_global_monthly_1950_2014.nc';

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
ncidSB = netcdf.create(file_tpb,'netcdf4');

lon_dim = netcdf.defDim(ncidSB,'lon',ni);
lat_dim = netcdf.defDim(ncidSB,'lat',nj);
time_dim = netcdf.defDim(ncidSB,'time',nt);

vidlat = netcdf.defVar(ncidSB,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSB,vidlat,'long_name','lat');
netcdf.putAtt(ncidSB,vidlat,'standard_name','lat');
netcdf.putAtt(ncidSB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidSB,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidSB,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSB,vidlon,'long_name','lon');
netcdf.putAtt(ncidSB,vidlon,'standard_name','lon');
netcdf.putAtt(ncidSB,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidSB,vidlon,'axis','X');

vidtSB = netcdf.defVar(ncidSB,'time','double',time_dim);
netcdf.putAtt(ncidSB,vidtSB,'long_name','time');
netcdf.putAtt(ncidSB,vidtSB,'standard_name','time');
netcdf.putAtt(ncidSB,vidtSB,'units','months since 1601-01-01' );
netcdf.putAtt(ncidSB,vidtSB,'calendar','360_day');
netcdf.putAtt(ncidSB,vidtSB,'axis','T');

vidbioSB = netcdf.defVar(ncidSB,'tpb','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSB,vidbioSB,'long_name','total pelagic biomass density');
netcdf.putAtt(ncidSB,vidbioSB,'units','grams wet weight m-2' );
netcdf.defVarFill(ncidSB,vidbioSB,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSB,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidSB,varid,'contact','C. Petrik');
netcdf.putAtt(ncidSB,varid,'institution','Texas A&M University');
netcdf.putAtt(ncidSB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidSB);

netcdf.putVar(ncidSB,vidlat,LAT);
netcdf.putVar(ncidSB,vidlon,LON);
netcdf.putVar(ncidSB,vidbioSB,tpb);
netcdf.putVar(ncidSB,vidtSB,htime);

netcdf.close(ncidSB);
%%
ncdisp(file_tpb)

%% tdb
ncidSD = netcdf.create(file_tdb,'netcdf4');

lon_dim = netcdf.defDim(ncidSD,'lon',ni);
lat_dim = netcdf.defDim(ncidSD,'lat',nj);
time_dim = netcdf.defDim(ncidSD,'time',nt);

vidlat = netcdf.defVar(ncidSD,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSD,vidlat,'long_name','lat');
netcdf.putAtt(ncidSD,vidlat,'standard_name','lat');
netcdf.putAtt(ncidSD,vidlat,'units','degrees_north');
netcdf.putAtt(ncidSD,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidSD,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSD,vidlon,'long_name','lon');
netcdf.putAtt(ncidSD,vidlon,'standard_name','lon');
netcdf.putAtt(ncidSD,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidSD,vidlon,'axis','X');

vidtSD = netcdf.defVar(ncidSD,'time','double',time_dim);
netcdf.putAtt(ncidSD,vidtSD,'long_name','time');
netcdf.putAtt(ncidSD,vidtSD,'standard_name','time');
netcdf.putAtt(ncidSD,vidtSD,'calendar','360_day');
netcdf.putAtt(ncidSD,vidtSD,'axis','T');
netcdf.putAtt(ncidSD,vidtSD,'units','months since 1601-01-01' );

vidbioSD = netcdf.defVar(ncidSD,'tdb','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSD,vidbioSD,'long_name','total demersal biomass density');
netcdf.putAtt(ncidSD,vidbioSD,'units','grams wet weight m-2' );
netcdf.defVarFill(ncidSD,vidbioSD,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSD,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSD,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidSD,varid,'contact','C. Petrik');
netcdf.putAtt(ncidSD,varid,'institution','Texas A&M University');
netcdf.putAtt(ncidSD,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidSD,varid,'includes benthos','no');

netcdf.endDef(ncidSD);

netcdf.putVar(ncidSD,vidlat,LAT);
netcdf.putVar(ncidSD,vidlon,LON);
netcdf.putVar(ncidSD,vidbioSD,tdb);
netcdf.putVar(ncidSD,vidtSD,htime);

netcdf.close(ncidSD);

%% tcb
ncidCB = netcdf.create(file_tcb,'netcdf4');

lon_dim = netcdf.defDim(ncidCB,'lon',ni);
lat_dim = netcdf.defDim(ncidCB,'lat',nj);
time_dim = netcdf.defDim(ncidCB,'time',nt);

vidlat = netcdf.defVar(ncidCB,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCB,vidlat,'long_name','lat');
netcdf.putAtt(ncidCB,vidlat,'standard_name','lat');
netcdf.putAtt(ncidCB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidCB,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidCB,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidCB,vidlon,'long_name','lon');
netcdf.putAtt(ncidCB,vidlon,'standard_name','lon');
netcdf.putAtt(ncidCB,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidCB,vidlon,'axis','X');

vidtCB = netcdf.defVar(ncidCB,'time','double',time_dim);
netcdf.putAtt(ncidCB,vidtCB,'long_name','time');
netcdf.putAtt(ncidCB,vidtCB,'standard_name','time');
netcdf.putAtt(ncidCB,vidtCB,'calendar','360_day');
netcdf.putAtt(ncidCB,vidtCB,'axis','T');
netcdf.putAtt(ncidCB,vidtCB,'units','months since 1601-01-01' );

vidbioCB = netcdf.defVar(ncidCB,'tcb','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidCB,vidbioCB,'long_name','total consumer biomass density');
netcdf.putAtt(ncidCB,vidbioCB,'units','grams wet weight m-2' );
netcdf.defVarFill(ncidCB,vidbioCB,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidCB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidCB,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidCB,varid,'contact','C. Petrik');
netcdf.putAtt(ncidCB,varid,'institution','Texas A&M University');
netcdf.putAtt(ncidCB,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidCB,varid,'includes benthos','yes');

netcdf.endDef(ncidCB);

netcdf.putVar(ncidCB,vidlat,LAT);
netcdf.putVar(ncidCB,vidlon,LON);
netcdf.putVar(ncidCB,vidbioCB,tcb);
netcdf.putVar(ncidCB,vidtCB,htime);

netcdf.close(ncidCB);

%% bp30cm
ncid30 = netcdf.create(file_bp30,'netcdf4');

lon_dim = netcdf.defDim(ncid30,'lon',ni);
lat_dim = netcdf.defDim(ncid30,'lat',nj);
time_dim = netcdf.defDim(ncid30,'time',nt);

vidlat = netcdf.defVar(ncid30,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid30,vidlat,'long_name','lat');
netcdf.putAtt(ncid30,vidlat,'standard_name','lat');
netcdf.putAtt(ncid30,vidlat,'units','degrees_north');
netcdf.putAtt(ncid30,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncid30,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid30,vidlon,'long_name','lon');
netcdf.putAtt(ncid30,vidlon,'standard_name','lon');
netcdf.putAtt(ncid30,vidlon,'units','degrees_east' );
netcdf.putAtt(ncid30,vidlon,'axis','X');

vidt30 = netcdf.defVar(ncid30,'time','double',time_dim);
netcdf.putAtt(ncid30,vidt30,'long_name','time');
netcdf.putAtt(ncid30,vidt30,'standard_name','time');
netcdf.putAtt(ncid30,vidt30,'calendar','360_day');
netcdf.putAtt(ncid30,vidt30,'axis','T');
netcdf.putAtt(ncid30,vidt30,'units','months since 1601-01-01' );

vidbio30 = netcdf.defVar(ncid30,'bp30cm','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncid30,vidbio30,'long_name','biomass density of pelagic < 30cm');
netcdf.putAtt(ncid30,vidbio30,'units','grams wet weight  m-2' );
netcdf.defVarFill(ncid30,vidbio30,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid30,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid30,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncid30,varid,'contact','C. Petrik');
netcdf.putAtt(ncid30,varid,'institution','Texas A&M University');
netcdf.putAtt(ncid30,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncid30);

netcdf.putVar(ncid30,vidlat,LAT);
netcdf.putVar(ncid30,vidlon,LON);
netcdf.putVar(ncid30,vidbio30,bp30cm);
netcdf.putVar(ncid30,vidt30,htime);

netcdf.close(ncid30);

%% bp90cm
ncid90 = netcdf.create(file_bp90,'netcdf4');

lon_dim = netcdf.defDim(ncid90,'lon',ni);
lat_dim = netcdf.defDim(ncid90,'lat',nj);
time_dim = netcdf.defDim(ncid90,'time',nt);

vidlat = netcdf.defVar(ncid90,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid90,vidlat,'long_name','lat');
netcdf.putAtt(ncid90,vidlat,'standard_name','lat');
netcdf.putAtt(ncid90,vidlat,'units','degrees_north');
netcdf.putAtt(ncid90,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncid90,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid90,vidlon,'long_name','lon');
netcdf.putAtt(ncid90,vidlon,'standard_name','lon');
netcdf.putAtt(ncid90,vidlon,'units','degrees_east' );
netcdf.putAtt(ncid90,vidlon,'axis','X');

vidt90 = netcdf.defVar(ncid90,'time','double',time_dim);
netcdf.putAtt(ncid90,vidt90,'long_name','time');
netcdf.putAtt(ncid90,vidt90,'standard_name','time');
netcdf.putAtt(ncid90,vidt90,'calendar','360_day');
netcdf.putAtt(ncid90,vidt90,'axis','T');
netcdf.putAtt(ncid90,vidt90,'units','months since 1601-01-01' );

vidbio90 = netcdf.defVar(ncid90,'bp90cm','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncid90,vidbio90,'long_name','biomass density of pelagic >=90cm');
netcdf.putAtt(ncid90,vidbio90,'units','grams wet weight  m-2' );
netcdf.defVarFill(ncid90,vidbio90,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid90,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid90,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncid90,varid,'contact','C. Petrik');
netcdf.putAtt(ncid90,varid,'institution','Texas A&M University');
netcdf.putAtt(ncid90,varid,'wet weight:C ratio','9:1');
% netcdf.putAtt(ncid90,varid,'feisty fish size','29.24 to 232.08cm');

netcdf.endDef(ncid90);

netcdf.putVar(ncid90,vidlat,LAT);
netcdf.putVar(ncid90,vidlon,LON);
netcdf.putVar(ncid90,vidbio90,bp90cm);
netcdf.putVar(ncid90,vidt90,htime);

netcdf.close(ncid90);

%% bd90cm
ncid90 = netcdf.create(file_bd90,'netcdf4');

lon_dim = netcdf.defDim(ncid90,'lon',ni);
lat_dim = netcdf.defDim(ncid90,'lat',nj);
time_dim = netcdf.defDim(ncid90,'time',nt);

vidlat = netcdf.defVar(ncid90,'lat','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid90,vidlat,'long_name','lat');
netcdf.putAtt(ncid90,vidlat,'standard_name','lat');
netcdf.putAtt(ncid90,vidlat,'units','degrees_north');
netcdf.putAtt(ncid90,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncid90,'lon','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncid90,vidlon,'long_name','lon');
netcdf.putAtt(ncid90,vidlon,'standard_name','lon');
netcdf.putAtt(ncid90,vidlon,'units','degrees_east' );
netcdf.putAtt(ncid90,vidlon,'axis','X');

vidt90 = netcdf.defVar(ncid90,'time','double',time_dim);
netcdf.putAtt(ncid90,vidt90,'long_name','time');
netcdf.putAtt(ncid90,vidt90,'standard_name','time');
netcdf.putAtt(ncid90,vidt90,'calendar','360_day');
netcdf.putAtt(ncid90,vidt90,'axis','T');
netcdf.putAtt(ncid90,vidt90,'units','months since 1601-01-01' );

vidbio90 = netcdf.defVar(ncid90,'bd90cm','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncid90,vidbio90,'long_name','biomass density of demersal >=90cm');
netcdf.putAtt(ncid90,vidbio90,'units','grams wet weight  m-2' );
netcdf.defVarFill(ncid90,vidbio90,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid90,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid90,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncid90,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncid90,varid,'contact','C. Petrik');
netcdf.putAtt(ncid90,varid,'institution','Texas A&M University');
netcdf.putAtt(ncid90,varid,'includes benthos','no');
% netcdf.putAtt(ncid90,varid,'feisty fish size','29.24 to 232.08cm');

netcdf.endDef(ncid90);

netcdf.putVar(ncid90,vidlat,LAT);
netcdf.putVar(ncid90,vidlon,LON);
netcdf.putVar(ncid90,vidbio90,bd90cm);
netcdf.putVar(ncid90,vidt90,htime);

netcdf.close(ncid90);


