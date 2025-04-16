% Calc Fish-MIP outputs saved as NetCDF
% Future time period
% SSP 534-over

clear 
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';
fpath=['/project/Feisty/NC/WG2300/' cfile '/IPSL/'];
harv = 'pristine';

%% Fish-MIP OUTPUTS =================================================

load([fpath 'FishMIP_outputs_monthly_IPSL_ssp534-over_pristine_' cfile '.mat'])

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
clear time

%% Grid
load('/project/Feisty/Fish-MIP/CMIP6/IPSL/Data_grid_ipsl_cmip6_2300.mat','GRD');
load('/project/Feisty/Fish-MIP/CMIP6/IPSL/gridspec_ipsl_cmip6_2300.mat');

[ni,nj] = size(LAT);

load('/project/Feisty/Fish-MIP/CMIP6/IPSL/ssp534over/ipsl_ssp534-over_time_monthly_2040_2300.mat');
t_all = time;

%% Reshape to lat,lon,yr
[nid,nt] = size(allC);

tpb = 1.000000020040877e20*ones(ni,nj,nt);
tdb = tpb;
tcb = tpb;
bp30cm = tpb;
bp90cm = tpb;
bd90cm = tpb;

for y=1:nt
    gtpb = 1.000000020040877e20*ones(ni,nj);
    ttpb = allPel(:,y);
    gtpb(GRD.ID) = ttpb;
    tpb(:,:,y) = gtpb;

    gtdb = 1.000000020040877e20*ones(ni,nj);
    ttdb = allD(:,y);
    gtdb(GRD.ID) = ttdb;
    tdb(:,:,y) = gtdb;

    gtcb = 1.000000020040877e20*ones(ni,nj);
    ttcb = allC(:,y);
    gtcb(GRD.ID) = ttcb;
    tcb(:,:,y) = gtcb;

    gp30cm = 1.000000020040877e20*ones(ni,nj);
    tp30cm = SPel(:,y);
    gp30cm(GRD.ID) = tp30cm;
    bp30cm(:,:,y) = gp30cm;

    gp90cm = 1.000000020040877e20*ones(ni,nj);
    tp90cm = LPel(:,y);
    gp90cm(GRD.ID) = tp90cm;
    bp90cm(:,:,y) = gp90cm;

    gd90cm = 1.000000020040877e20*ones(ni,nj);
    td90cm = allD(:,y);
    gd90cm(GRD.ID) = td90cm;
    bd90cm(:,:,y) = gd90cm;
end

%% Quick look
pb = tpb(:,:,105);
db = tdb(:,:,105);
cb = tcb(:,:,105);
bp30 = bp30cm(:,:,105);
bp90 = bp90cm(:,:,105);
bd90 = bd90cm(:,:,105);


%% is lat & lon correct?
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% Output data naming conventions
%<model>_<climate-forcing>_<bias-adjustment>_<climate-scenario>_
%<soc-scenario>_<sens-scenario>_<variable>_<global>_<timestep>_<start-year>_
%<end-year>.nc

%e.g.
%apecosm_ipsl-esm4_nobasd_picontrol_histsoc_default_tcb_global_monthly_2001_2010.nc

%% Setup netcdf path to store to
fname1 = 'feisty_ipsl-cm6a-lr_nobasd_ssp534-over_nat_default_';
fname2 = '_global_monthly_2040_2300.nc';

file_tpb = [fpath fname1 'tpb' fname2];
file_tdb = [fpath fname1 'tdb' fname2];
file_tcb = [fpath fname1 'tcb' fname2];
file_bp30 = [fpath fname1 'bp30cm' fname2];
file_bp90 = [fpath fname1 'bp90cm' fname2];
file_bd90 = [fpath fname1 'bd90cm' fname2];

[ni,nj,nt] = size(tpb);

%% Lat & Lon should be vectors
LAT = (LAT(1,:));
LON = LON(:,1);

%tpb.dtype="float64" should be "float32"
%'NC_FLOAT'

%tpb.chunking=contiguous should be [1, 180, 360] (with proper depencency order)

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
netcdf.putAtt(ncidSB,vidtSB,'units','months since 1601-01-01' );
netcdf.putAtt(ncidSB,vidtSB,'calendar','360_day');
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
netcdf.defVarChunking(ncidSB,vidbioSB,'CHUNKED',[1, 180, 360]);
netcdf.putAtt(ncidSB,vidbioSB,'long_name','Total Pelagic Biomass Density');
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
netcdf.putAtt(ncidSD,vidtSD,'calendar','360_day');
netcdf.putAtt(ncidSD,vidtSD,'axis','T');
netcdf.putAtt(ncidSD,vidtSD,'units','months since 1601-01-01' );

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
netcdf.putAtt(ncidSD,vidbioSD,'long_name','Total Demersal Biomass Density');
netcdf.putAtt(ncidSD,vidbioSD,'units','g m-2' );
netcdf.defVarFill(ncidSD,vidbioSD,false,1.000000020040877e+20);
netcdf.putAtt(ncidSD,vidbioSD,'missing value',1.000000020040877e+20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSD,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSD,varid,'_FillValue',1.000000020040877e+20);
netcdf.putAtt(ncidSB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidSB,varid,'institution','UCSD');
netcdf.putAtt(ncidSD,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidSD,varid,'includes benthos','no');

netcdf.endDef(ncidSD);

tdb = single(tdb);
netcdf.putVar(ncidSD,vidlat,LAT);
netcdf.putVar(ncidSD,vidlon,LON);
netcdf.putVar(ncidSD,vidbioSD,tdb);
netcdf.putVar(ncidSD,vidtSD,time);

netcdf.close(ncidSD);

%% tcb
ncidCB = netcdf.create(file_tcb,cmode);

time_dim = netcdf.defDim(ncidCB,'time',nt);
lon_dim = netcdf.defDim(ncidCB,'lon',ni);
lat_dim = netcdf.defDim(ncidCB,'lat',nj);

vidtCB = netcdf.defVar(ncidCB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidCB,vidtCB,'long_name','time');
netcdf.putAtt(ncidCB,vidtCB,'standard_name','time');
netcdf.putAtt(ncidCB,vidtCB,'calendar','360_day');
netcdf.putAtt(ncidCB,vidtCB,'axis','T');
netcdf.putAtt(ncidCB,vidtCB,'units','months since 1601-01-01' );

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
netcdf.putAtt(ncidCB,vidbioCB,'long_name','Total Consumer Biomass Density');
netcdf.putAtt(ncidCB,vidbioCB,'units','g m-2' );
netcdf.defVarFill(ncidCB,vidbioCB,false,1.000000020040877e+20);
netcdf.putAtt(ncidCB,vidbioCB,'missing value',1.000000020040877e+20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidCB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidCB,varid,'_FillValue',1.000000020040877e+20);
netcdf.putAtt(ncidSB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidSB,varid,'institution','UCSD');
netcdf.putAtt(ncidCB,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidCB,varid,'includes benthos','yes');

netcdf.endDef(ncidCB);

tcb = single(tcb);
netcdf.putVar(ncidCB,vidlat,LAT);
netcdf.putVar(ncidCB,vidlon,LON);
netcdf.putVar(ncidCB,vidbioCB,tcb);
netcdf.putVar(ncidCB,vidtCB,time);

netcdf.close(ncidCB);

%% bp30cm
ncid30 = netcdf.create(file_bp30,cmode);

time_dim = netcdf.defDim(ncid30,'time',nt);
lon_dim = netcdf.defDim(ncid30,'lon',ni);
lat_dim = netcdf.defDim(ncid30,'lat',nj);

vidt30 = netcdf.defVar(ncid30,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncid30,vidt30,'long_name','time');
netcdf.putAtt(ncid30,vidt30,'standard_name','time');
netcdf.putAtt(ncid30,vidt30,'calendar','360_day');
netcdf.putAtt(ncid30,vidt30,'axis','T');
netcdf.putAtt(ncid30,vidt30,'units','months since 1601-01-01' );

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
netcdf.putAtt(ncid30,vidbio30,'long_name','Biomass Density of Small Pelagics < 30cm');
netcdf.putAtt(ncid30,vidbio30,'units','g  m-2' );
netcdf.defVarFill(ncid30,vidbio30,false,1.000000020040877e+20);
netcdf.putAtt(ncid30,vidbio30,'missing value',1.000000020040877e+20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid30,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid30,varid,'_FillValue',1.000000020040877e+20);
netcdf.putAtt(ncidSB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidSB,varid,'institution','UCSD');
netcdf.putAtt(ncid30,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncid30);

bp30cm = single(bp30cm);
netcdf.putVar(ncid30,vidlat,LAT);
netcdf.putVar(ncid30,vidlon,LON);
netcdf.putVar(ncid30,vidbio30,bp30cm);
netcdf.putVar(ncid30,vidt30,time);

netcdf.close(ncid30);

%% bp90cm
ncid90 = netcdf.create(file_bp90,cmode);

time_dim = netcdf.defDim(ncid90,'time',nt);
lon_dim = netcdf.defDim(ncid90,'lon',ni);
lat_dim = netcdf.defDim(ncid90,'lat',nj);

vidt90 = netcdf.defVar(ncid90,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncid90,vidt90,'long_name','time');
netcdf.putAtt(ncid90,vidt90,'standard_name','time');
netcdf.putAtt(ncid90,vidt90,'calendar','360_day');
netcdf.putAtt(ncid90,vidt90,'axis','T');
netcdf.putAtt(ncid90,vidt90,'units','months since 1601-01-01' );

vidlon = netcdf.defVar(ncid90,'lon','NC_DOUBLE',lon_dim);
netcdf.putAtt(ncid90,vidlon,'long_name','longitude');
netcdf.putAtt(ncid90,vidlon,'standard_name','longitude');
netcdf.putAtt(ncid90,vidlon,'units','degrees_east' );
netcdf.putAtt(ncid90,vidlon,'axis','X');

vidlat = netcdf.defVar(ncid90,'lat','NC_DOUBLE',lat_dim);
netcdf.putAtt(ncid90,vidlat,'long_name','latitude');
netcdf.putAtt(ncid90,vidlat,'standard_name','latitude');
netcdf.putAtt(ncid90,vidlat,'units','degrees_north');
netcdf.putAtt(ncid90,vidlat,'axis','Y');

vidbio90 = netcdf.defVar(ncid90,'bp90cm','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncid90,vidbio90,'long_name','Biomass Density of Large Pelagics >=90cm');
netcdf.putAtt(ncid90,vidbio90,'units','g  m-2' );
netcdf.defVarFill(ncid90,vidbio90,false,1.000000020040877e+20);
netcdf.putAtt(ncid90,vidbio90,'missing value',1.000000020040877e+20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid90,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid90,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidSB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidSB,varid,'institution','UCSD');
netcdf.putAtt(ncid90,varid,'wet weight:C ratio','9:1');
% netcdf.putAtt(ncid90,varid,'feisty fish size','29.24 to 232.08cm');

netcdf.endDef(ncid90);

bp90cm = single(bp90cm);
netcdf.putVar(ncid90,vidlat,LAT);
netcdf.putVar(ncid90,vidlon,LON);
netcdf.putVar(ncid90,vidbio90,bp90cm);
netcdf.putVar(ncid90,vidt90,time);

netcdf.close(ncid90);

%% bd90cm
ncid90 = netcdf.create(file_bd90,cmode);

time_dim = netcdf.defDim(ncid90,'time',nt);
lon_dim = netcdf.defDim(ncid90,'lon',ni);
lat_dim = netcdf.defDim(ncid90,'lat',nj);

vidt90 = netcdf.defVar(ncid90,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncid90,vidt90,'long_name','time');
netcdf.putAtt(ncid90,vidt90,'standard_name','time');
netcdf.putAtt(ncid90,vidt90,'calendar','360_day');
netcdf.putAtt(ncid90,vidt90,'axis','T');
netcdf.putAtt(ncid90,vidt90,'units','months since 1601-01-01' );

vidlon = netcdf.defVar(ncid90,'lon','NC_DOUBLE',lon_dim);
netcdf.putAtt(ncid90,vidlon,'long_name','longitude');
netcdf.putAtt(ncid90,vidlon,'standard_name','longitude');
netcdf.putAtt(ncid90,vidlon,'units','degrees_east' );
netcdf.putAtt(ncid90,vidlon,'axis','X');

vidlat = netcdf.defVar(ncid90,'lat','NC_DOUBLE',lat_dim);
netcdf.putAtt(ncid90,vidlat,'long_name','latitude');
netcdf.putAtt(ncid90,vidlat,'standard_name','latitude');
netcdf.putAtt(ncid90,vidlat,'units','degrees_north');
netcdf.putAtt(ncid90,vidlat,'axis','Y');

vidbio90 = netcdf.defVar(ncid90,'bd90cm','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncid90,vidbio90,'long_name','Biomass Density of Large Demersals >=90cm');
netcdf.putAtt(ncid90,vidbio90,'units','g  m-2' );
netcdf.defVarFill(ncid90,vidbio90,false,1.000000020040877e+20);
netcdf.putAtt(ncid90,vidbio90,'missing value',1.000000020040877e+20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid90,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid90,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncid90,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidSB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidSB,varid,'institution','UCSD');
netcdf.putAtt(ncid90,varid,'includes benthos','no');
% netcdf.putAtt(ncid90,varid,'feisty fish size','29.24 to 232.08cm');

netcdf.endDef(ncid90);

bd90cm = single(bd90cm);
netcdf.putVar(ncid90,vidlat,LAT);
netcdf.putVar(ncid90,vidlon,LON);
netcdf.putVar(ncid90,vidbio90,bd90cm);
netcdf.putVar(ncid90,vidt90,time);

netcdf.close(ncid90);
