% FEISTY output at all locations
% Hist ctrlclim fished 1/4 degree
% Early periods (1841-1960) can have 1841 as the reference year 
% Historic experimental period uses 1901 as the reference year

clear 
close all

load('FishMIP_phase3a_exper_times.mat')
time_long_name = 'time';
time_standard_name = 'time';
time_units = 'days since 1901-1-1 00:00:00'; %hist_units;
time_axis = 'T';
calendar = '365_day';

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];

%% SD
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_v3.2_empHP_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[nid,nt] = size(biomass);

SD.bio = biomass;
clear biomass 

% MD
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_v3.2_empHP_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
MD.yield = yield;
clear yield
clear biomass

% LD
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_v3.2_empHP_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
LD.yield = yield;
clear yield
clear biomass

allDB = SD.bio + MD.bio + LD.bio;
allDC = MD.yield + LD.yield;

clear SD MD LD

%% catch totals per month
% mult mean catch per day by # days each mo
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
mos = repmat(MNTH,nid,(nt/12));

allDC = allDC .* mos;

%% ========================== netcdf info ===============================

%% ESM
epath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
load([epath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'],'GRD');
load([epath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'],'LAT','LON');

[ni,nj] = size(LAT);

lat = LAT(1,:);
lon = LON(:,1);

%% Reshape to lat,lon,yr
[nid,nt] = size(allDB);

tpb = 1.0e20*ones(ni,nj,nt);
tdb = tpb;
tdc = tpb;

for y=1:nt
    gtdb = 1.0e20*ones(ni,nj);
    ttdb = allDB(:,y);
    gtdb(GRD.ID) = ttdb;
    tdb(:,:,y) = gtdb;
    
    gtdc = 1.0e20*ones(ni,nj);
    ttdc = allDC(:,y);
    gtdc(GRD.ID) = ttdc;
    tdc(:,:,y) = gtdc;
    
end

%%
clear allDB allDC 

%% Quick look
db = tdb(:,:,150);

figure(2)
pcolor(log10(db))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('all Dem')

%% Output data naming conventions
%<model>_<climate-forcing>_<climate-scenario>_<soc-scenario>_
% <sens-scenario>_<variable>_<region>_<time-step>_<start-year>_
% <end-year>.nc

% climate scenario:   obsclim or ctrlclim
% socioecon scenario: histsoc or nat
% sens scenario:      15arcmin or onedeg

%e.g.
%boats_gfdl-mom6_cobalt2_none_obsclim_histsoc_default_tcb_global_monthly_1840_2010.nc

close all

%% Setup netcdf path to store to
fname1 = 'feisty_gfdl-mom6-cobalt2_obsclim_histsoc_1955-riverine-input_';
fname2 = '_global_monthly_1961_2010.nc';

file_tdb = [fpath fname1 'tdb' fname2];
file_bd90 = [fpath fname1 'bd90cm' fname2];

file_tdc = [fpath fname1 'tdc' fname2];
file_cd90 = [fpath fname1 'cd90cm' fname2];

[ni,nj,nt] = size(tpb);

cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

%% tdb
ncidSD = netcdf.create(file_tdb,cmode);

lon_dim = netcdf.defDim(ncidSD,'lon',ni);
lat_dim = netcdf.defDim(ncidSD,'lat',nj);
time_dim = netcdf.defDim(ncidSD,'time',nt);

vidlat = netcdf.defVar(ncidSD,'lat','double',lat_dim);
netcdf.putAtt(ncidSD,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSD,vidlat,'standard_name','lat');
netcdf.putAtt(ncidSD,vidlat,'units','degrees_north');
netcdf.putAtt(ncidSD,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidSD,'lon','double',lon_dim);
netcdf.putAtt(ncidSD,vidlon,'long_name','longitude');
netcdf.putAtt(ncidSD,vidlon,'standard_name','lon');
netcdf.putAtt(ncidSD,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidSD,vidlon,'axis','X');

vidtSD = netcdf.defVar(ncidSD,'time','double',time_dim);
netcdf.putAtt(ncidSD,vidtSD,'long_name',time_long_name);
netcdf.putAtt(ncidSD,vidtSD,'long_name',time_standard_name);
netcdf.putAtt(ncidSD,vidtSD,'calendar',calendar);
netcdf.putAtt(ncidSD,vidtSD,'axis',time_axis);
netcdf.putAtt(ncidSD,vidtSD,'units',time_units);

vidbioSD = netcdf.defVar(ncidSD,'tdb','double',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidSD,vidbioSD,'CHUNKED',[ni nj 1]);
netcdf.putAtt(ncidSD,vidbioSD,'long_name','Total Demersal Biomass Density');
netcdf.putAtt(ncidSD,vidbioSD,'units','g m-2' );
netcdf.defVarFill(ncidSD,vidbioSD,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSD,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSD,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidSD,varid,'contact','C. Petrik');
netcdf.putAtt(ncidSD,varid,'institution','UC San Diego');
netcdf.putAtt(ncidSD,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidSD,varid,'includes benthos','no');

netcdf.endDef(ncidSD);

netcdf.putVar(ncidSD,vidlat,lat);
netcdf.putVar(ncidSD,vidlon,lon);
netcdf.putVar(ncidSD,vidbioSD,tdb);
netcdf.putVar(ncidSD,vidtSD,hist_time);

netcdf.close(ncidSD);

%% bd90cm
ncid90 = netcdf.create(file_bd90,cmode);

lon_dim = netcdf.defDim(ncid90,'lon',ni);
lat_dim = netcdf.defDim(ncid90,'lat',nj);
time_dim = netcdf.defDim(ncid90,'time',nt);

vidlat = netcdf.defVar(ncid90,'lat','double',lat_dim);
netcdf.putAtt(ncid90,vidlat,'long_name','latitude');
netcdf.putAtt(ncid90,vidlat,'standard_name','lat');
netcdf.putAtt(ncid90,vidlat,'units','degrees_north');
netcdf.putAtt(ncid90,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncid90,'lon','double',lon_dim);
netcdf.putAtt(ncid90,vidlon,'long_name','longitude');
netcdf.putAtt(ncid90,vidlon,'standard_name','lon');
netcdf.putAtt(ncid90,vidlon,'units','degrees_east' );
netcdf.putAtt(ncid90,vidlon,'axis','X');

vidt90 = netcdf.defVar(ncid90,'time','double',time_dim);
netcdf.putAtt(ncid90,vidt90,'long_name',time_long_name);
netcdf.putAtt(ncid90,vidt90,'long_name',time_standard_name);
netcdf.putAtt(ncid90,vidt90,'calendar',calendar);
netcdf.putAtt(ncid90,vidt90,'axis',time_axis);
netcdf.putAtt(ncid90,vidt90,'units',time_units);

vidbio90 = netcdf.defVar(ncid90,'bd90cm','double',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncid90,vidbio90,'CHUNKED',[ni nj 1]);
netcdf.putAtt(ncid90,vidbio90,'long_name','Biomass Density of Large Demersals >=90cm');
netcdf.putAtt(ncid90,vidbio90,'units','g m-2' );
netcdf.defVarFill(ncid90,vidbio90,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid90,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid90,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncid90,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncid90,varid,'contact','C. Petrik');
netcdf.putAtt(ncid90,varid,'institution','UC San Diego');
netcdf.putAtt(ncid90,varid,'includes benthos','no');
% netcdf.putAtt(ncid90,varid,'feisty fish size','29.24 to 232.08cm');

netcdf.endDef(ncid90);

netcdf.putVar(ncid90,vidlat,lat);
netcdf.putVar(ncid90,vidlon,lon);
netcdf.putVar(ncid90,vidbio90,tdb);
netcdf.putVar(ncid90,vidt90,hist_time);

netcdf.close(ncid90);


%%  -------------------------CATCH-----------------------------


%% tdc
ncidSD = netcdf.create(file_tdc,cmode);

lon_dim = netcdf.defDim(ncidSD,'lon',ni);
lat_dim = netcdf.defDim(ncidSD,'lat',nj);
time_dim = netcdf.defDim(ncidSD,'time',nt);

vidlat = netcdf.defVar(ncidSD,'lat','double',lat_dim);
netcdf.putAtt(ncidSD,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSD,vidlat,'standard_name','lat');
netcdf.putAtt(ncidSD,vidlat,'units','degrees_north');
netcdf.putAtt(ncidSD,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidSD,'lon','double',lon_dim);
netcdf.putAtt(ncidSD,vidlon,'long_name','longitude');
netcdf.putAtt(ncidSD,vidlon,'standard_name','lon');
netcdf.putAtt(ncidSD,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidSD,vidlon,'axis','X');

vidtSD = netcdf.defVar(ncidSD,'time','double',time_dim);
netcdf.putAtt(ncidSD,vidtSD,'long_name',time_long_name);
netcdf.putAtt(ncidSD,vidtSD,'long_name',time_standard_name);
netcdf.putAtt(ncidSD,vidtSD,'calendar',calendar);
netcdf.putAtt(ncidSD,vidtSD,'axis',time_axis);
netcdf.putAtt(ncidSD,vidtSD,'units',time_units);

vidbioSD = netcdf.defVar(ncidSD,'tdc','double',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidSD,vidbioSD,'CHUNKED',[ni nj 1]);
netcdf.putAtt(ncidSD,vidbioSD,'long_name','Total Demersal Catch');
netcdf.putAtt(ncidSD,vidbioSD,'units','g m-2' );
netcdf.defVarFill(ncidSD,vidbioSD,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSD,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSD,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidSD,varid,'contact','C. Petrik');
netcdf.putAtt(ncidSD,varid,'institution','UC San Diego');
netcdf.putAtt(ncidSD,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidSD,varid,'includes benthos','no');

netcdf.endDef(ncidSD);

netcdf.putVar(ncidSD,vidlat,lat);
netcdf.putVar(ncidSD,vidlon,lon);
netcdf.putVar(ncidSD,vidbioSD,tdc);
netcdf.putVar(ncidSD,vidtSD,hist_time);

netcdf.close(ncidSD);

%% cd90cm
ncid90 = netcdf.create(file_cd90,cmode);

lon_dim = netcdf.defDim(ncid90,'lon',ni);
lat_dim = netcdf.defDim(ncid90,'lat',nj);
time_dim = netcdf.defDim(ncid90,'time',nt);

vidlat = netcdf.defVar(ncid90,'lat','double',lat_dim);
netcdf.putAtt(ncid90,vidlat,'long_name','latitude');
netcdf.putAtt(ncid90,vidlat,'standard_name','lat');
netcdf.putAtt(ncid90,vidlat,'units','degrees_north');
netcdf.putAtt(ncid90,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncid90,'lon','double',lon_dim);
netcdf.putAtt(ncid90,vidlon,'long_name','longitude');
netcdf.putAtt(ncid90,vidlon,'standard_name','lon');
netcdf.putAtt(ncid90,vidlon,'units','degrees_east' );
netcdf.putAtt(ncid90,vidlon,'axis','X');

vidt90 = netcdf.defVar(ncid90,'time','double',time_dim);
netcdf.putAtt(ncid90,vidt90,'long_name',time_long_name);
netcdf.putAtt(ncid90,vidt90,'long_name',time_standard_name);
netcdf.putAtt(ncid90,vidt90,'calendar',calendar);
netcdf.putAtt(ncid90,vidt90,'axis',time_axis);
netcdf.putAtt(ncid90,vidt90,'units',time_units);

vidbio90 = netcdf.defVar(ncid90,'cd90cm','double',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncid90,vidbio90,'CHUNKED',[ni nj 1]);
netcdf.putAtt(ncid90,vidbio90,'long_name','Catch Density of Large Demersals >=90cm');
netcdf.putAtt(ncid90,vidbio90,'units','g m-2' );
netcdf.defVarFill(ncid90,vidbio90,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid90,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid90,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncid90,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncid90,varid,'contact','C. Petrik');
netcdf.putAtt(ncid90,varid,'institution','UC San Diego');
netcdf.putAtt(ncid90,varid,'includes benthos','no');
% netcdf.putAtt(ncid90,varid,'feisty fish size','29.24 to 232.08cm');

netcdf.endDef(ncid90);

netcdf.putVar(ncid90,vidlat,lat);
netcdf.putVar(ncid90,vidlon,lon);
netcdf.putVar(ncid90,vidbio90,tdc);
netcdf.putVar(ncid90,vidt90,hist_time);

netcdf.close(ncid90);







