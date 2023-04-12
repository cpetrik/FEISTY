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

%% SP
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_v3.2_empHP_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);
[nid,nt] = size(biomass);

SP.bio = biomass;
clear biomass

% MP
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_v3.2_empHP_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
MP.yield = yield;
clear yield
clear biomass

% LP
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_v3.2_empHP_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
LP.yield = yield;
clear yield
clear biomass

allPB = SP.bio + MP.bio + LP.bio;
allPC = MP.yield + LP.yield;

clear SP MP LP

%% SF 
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_v3.2_empHP_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
clear biomass 

% MF
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_v3.2_empHP_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
MF.yield = yield;
clear yield
clear biomass

allFB = SF.bio + MF.bio;
allFC = MF.yield;

clear SF MF

%% catch totals per month
% mult mean catch per day by # days each mo
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
mos = repmat(MNTH,nid,(nt/12));

allFC = allFC .* mos;
allPC = allPC .* mos;

%% ========================== netcdf info ===============================

%% ESM
epath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
load([epath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'],'GRD');
load([epath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'],'LAT','LON');

[ni,nj] = size(LAT);

lat = LAT(1,:);
lon = LON(:,1);

%% Reshape to lat,lon,yr
[nid,nt] = size(allPB);

tpb = 1.0e20*ones(ni,nj,nt);
tfb = tpb;
tfc = tpb;
tpc = tpb;

for y=1:nt
    gtpb = 1.0e20*ones(ni,nj);
    ttpb = allPB(:,y);
    gtpb(GRD.ID) = ttpb;
    tpb(:,:,y) = gtpb;
    
    gtfb = 1.0e20*ones(ni,nj);
    ttfb = allFB(:,y);
    gtfb(GRD.ID) = ttfb;
    tfb(:,:,y) = gtfb;

    gtpc = 1.0e20*ones(ni,nj);
    ttpc = allPC(:,y);
    gtpc(GRD.ID) = ttpc;
    tpc(:,:,y) = gtpc;
    
    gtfc = 1.0e20*ones(ni,nj);
    ttfc = allFC(:,y);
    gtfc(GRD.ID) = ttfc;
    tfc(:,:,y) = gtfc;
end

%%
clear allPB allFB allPC allFC

%%
allPelB = tfb + tpb;

allPelC = tfc + tpc;

%% Quick look
pb = tpb(:,:,150);
fb = tfb(:,:,150);

figure(1)
pcolor(log10(pb))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('all LPel')

figure(3)
pcolor(log10(fb))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('all F')

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

file_tpb = [fpath fname1 'tpb' fname2];
file_bp30 = [fpath fname1 'bp30cm' fname2];
file_bp90 = [fpath fname1 'bp90cm' fname2];

file_tpc = [fpath fname1 'tpc' fname2];
file_cp30 = [fpath fname1 'cp30cm' fname2];
file_cp90 = [fpath fname1 'cp90cm' fname2];

[ni,nj,nt] = size(tpb);

cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

%% tpb
%ncidSB = netcdf.create(file_tpb,'netcdf4');
ncidSB = netcdf.create(file_tpb,cmode);

lon_dim = netcdf.defDim(ncidSB,'lon',ni);
lat_dim = netcdf.defDim(ncidSB,'lat',nj);
time_dim = netcdf.defDim(ncidSB,'time',nt);

vidlat = netcdf.defVar(ncidSB,'lat','double',lat_dim);
netcdf.putAtt(ncidSB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'standard_name','lat');
netcdf.putAtt(ncidSB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidSB,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidSB,'lon','double',lon_dim);
netcdf.putAtt(ncidSB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidSB,vidlon,'standard_name','lon');
netcdf.putAtt(ncidSB,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidSB,vidlon,'axis','X');

vidtSB = netcdf.defVar(ncidSB,'time','double',time_dim);
netcdf.putAtt(ncidSB,vidtSB,'long_name',time_long_name);
netcdf.putAtt(ncidSB,vidtSB,'standard_name',time_standard_name);
netcdf.putAtt(ncidSB,vidtSB,'units',time_units);
netcdf.putAtt(ncidSB,vidtSB,'calendar',calendar);
netcdf.putAtt(ncidSB,vidtSB,'axis',time_axis);

vidbioSB = netcdf.defVar(ncidSB,'tpb','double',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidSB,vidbioSB,'CHUNKED',[ni nj 1]);
netcdf.putAtt(ncidSB,vidbioSB,'long_name','Total Pelagic Biomass Density');
netcdf.putAtt(ncidSB,vidbioSB,'units','g m-2' );
netcdf.defVarFill(ncidSB,vidbioSB,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSB,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidSB,varid,'contact','C. Petrik');
netcdf.putAtt(ncidSB,varid,'institution','UC San Diego');
netcdf.putAtt(ncidSB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidSB);

netcdf.putVar(ncidSB,vidlat,lat);
netcdf.putVar(ncidSB,vidlon,lon);
netcdf.putVar(ncidSB,vidbioSB,allPelB);
netcdf.putVar(ncidSB,vidtSB,hist_time);

netcdf.close(ncidSB);
%%
ncdisp(file_tpb)

%% bp30cm
ncid30 = netcdf.create(file_bp30,cmode);

lon_dim = netcdf.defDim(ncid30,'lon',ni);
lat_dim = netcdf.defDim(ncid30,'lat',nj);
time_dim = netcdf.defDim(ncid30,'time',nt);

vidlat = netcdf.defVar(ncid30,'lat','double',lat_dim);
netcdf.putAtt(ncid30,vidlat,'long_name','latitude');
netcdf.putAtt(ncid30,vidlat,'standard_name','lat');
netcdf.putAtt(ncid30,vidlat,'units','degrees_north');
netcdf.putAtt(ncid30,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncid30,'lon','double',lon_dim);
netcdf.putAtt(ncid30,vidlon,'long_name','longitude');
netcdf.putAtt(ncid30,vidlon,'standard_name','lon');
netcdf.putAtt(ncid30,vidlon,'units','degrees_east' );
netcdf.putAtt(ncid30,vidlon,'axis','X');

vidt30 = netcdf.defVar(ncid30,'time','double',time_dim);
netcdf.putAtt(ncid30,vidt30,'long_name',time_long_name);
netcdf.putAtt(ncid30,vidt30,'long_name',time_standard_name);
netcdf.putAtt(ncid30,vidt30,'calendar',calendar);
netcdf.putAtt(ncid30,vidt30,'axis',time_axis);
netcdf.putAtt(ncid30,vidt30,'units',time_units);

vidbio30 = netcdf.defVar(ncid30,'bp30cm','double',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncid30,vidbio30,'CHUNKED',[ni nj 1]);
netcdf.putAtt(ncid30,vidbio30,'long_name','Biomass Density of Small Pelagics <30cm');
netcdf.putAtt(ncid30,vidbio30,'units','g m-2' );
netcdf.defVarFill(ncid30,vidbio30,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid30,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid30,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncid30,varid,'contact','C. Petrik');
netcdf.putAtt(ncid30,varid,'institution','UC San Diego');
netcdf.putAtt(ncid30,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncid30);

netcdf.putVar(ncid30,vidlat,lat);
netcdf.putVar(ncid30,vidlon,lon);
netcdf.putVar(ncid30,vidbio30,tfb);
netcdf.putVar(ncid30,vidt30,hist_time);

netcdf.close(ncid30);

%% bp90cm
ncid90 = netcdf.create(file_bp90,cmode);

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

vidbio90 = netcdf.defVar(ncid90,'bp90cm','double',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncid90,vidbio90,'CHUNKED',[ni nj 1]);
netcdf.putAtt(ncid90,vidbio90,'long_name','Biomass Density of Large Pelagics >=90cm');
netcdf.putAtt(ncid90,vidbio90,'units','g m-2' );
netcdf.defVarFill(ncid90,vidbio90,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid90,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid90,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncid90,varid,'contact','C. Petrik');
netcdf.putAtt(ncid90,varid,'institution','UC San Diego');
netcdf.putAtt(ncid90,varid,'wet weight:C ratio','9:1');
% netcdf.putAtt(ncid90,varid,'feisty fish size','29.24 to 232.08cm');

netcdf.endDef(ncid90);

netcdf.putVar(ncid90,vidlat,lat);
netcdf.putVar(ncid90,vidlon,lon);
netcdf.putVar(ncid90,vidbio90,tpb);
netcdf.putVar(ncid90,vidt90,hist_time);

netcdf.close(ncid90);

%% tpc -------------------------CATCH-----------------------------
ncidSB = netcdf.create(file_tpc,cmode);

lon_dim = netcdf.defDim(ncidSB,'lon',ni);
lat_dim = netcdf.defDim(ncidSB,'lat',nj);
time_dim = netcdf.defDim(ncidSB,'time',nt);

vidlat = netcdf.defVar(ncidSB,'lat','double',lat_dim);
netcdf.putAtt(ncidSB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'standard_name','lat');
netcdf.putAtt(ncidSB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidSB,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidSB,'lon','double',lon_dim);
netcdf.putAtt(ncidSB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidSB,vidlon,'standard_name','lon');
netcdf.putAtt(ncidSB,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidSB,vidlon,'axis','X');

vidtSB = netcdf.defVar(ncidSB,'time','double',time_dim);
netcdf.putAtt(ncidSB,vidtSB,'long_name',time_long_name);
netcdf.putAtt(ncidSB,vidtSB,'long_name',time_standard_name);
netcdf.putAtt(ncidSB,vidtSB,'units',time_units);
netcdf.putAtt(ncidSB,vidtSB,'calendar',calendar);
netcdf.putAtt(ncidSB,vidtSB,'axis',time_axis);

vidbioSB = netcdf.defVar(ncidSB,'tpc','double',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidSB,vidbioSB,'CHUNKED',[ni nj 1]);
netcdf.putAtt(ncidSB,vidbioSB,'long_name','Total Pelagic Catch');
netcdf.putAtt(ncidSB,vidbioSB,'units','g m-2' );
netcdf.defVarFill(ncidSB,vidbioSB,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSB,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidSB,varid,'contact','C. Petrik');
netcdf.putAtt(ncidSB,varid,'institution','UC San Diego');
netcdf.putAtt(ncidSB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidSB);

netcdf.putVar(ncidSB,vidlat,lat);
netcdf.putVar(ncidSB,vidlon,lon);
netcdf.putVar(ncidSB,vidbioSB,allPelC);
netcdf.putVar(ncidSB,vidtSB,hist_time);

netcdf.close(ncidSB);

%% cp30cm
ncid30 = netcdf.create(file_cp30,cmode);

lon_dim = netcdf.defDim(ncid30,'lon',ni);
lat_dim = netcdf.defDim(ncid30,'lat',nj);
time_dim = netcdf.defDim(ncid30,'time',nt);

vidlat = netcdf.defVar(ncid30,'lat','double',lat_dim);
netcdf.putAtt(ncid30,vidlat,'long_name','latitude');
netcdf.putAtt(ncid30,vidlat,'standard_name','lat');
netcdf.putAtt(ncid30,vidlat,'units','degrees_north');
netcdf.putAtt(ncid30,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncid30,'lon','double',lon_dim);
netcdf.putAtt(ncid30,vidlon,'long_name','longitude');
netcdf.putAtt(ncid30,vidlon,'standard_name','lon');
netcdf.putAtt(ncid30,vidlon,'units','degrees_east' );
netcdf.putAtt(ncid30,vidlon,'axis','X');

vidt30 = netcdf.defVar(ncid30,'time','double',time_dim);
netcdf.putAtt(ncid30,vidt30,'long_name',time_long_name);
netcdf.putAtt(ncid30,vidt30,'long_name',time_standard_name);
netcdf.putAtt(ncid30,vidt30,'calendar',calendar);
netcdf.putAtt(ncid30,vidt30,'axis',time_axis);
netcdf.putAtt(ncid30,vidt30,'units',time_units);

vidbio30 = netcdf.defVar(ncid30,'cp30cm','double',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncid30,vidbio30,'CHUNKED',[ni nj 1]);
netcdf.putAtt(ncid30,vidbio30,'long_name','Catch Density of Small Pelagics <30cm');
netcdf.putAtt(ncid30,vidbio30,'units','g m-2' );
netcdf.defVarFill(ncid30,vidbio30,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid30,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid30,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncid30,varid,'contact','C. Petrik');
netcdf.putAtt(ncid30,varid,'institution','UC San Diego');
netcdf.putAtt(ncid30,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncid30);

netcdf.putVar(ncid30,vidlat,lat);
netcdf.putVar(ncid30,vidlon,lon);
netcdf.putVar(ncid30,vidbio30,tfc);
netcdf.putVar(ncid30,vidt30,hist_time);

netcdf.close(ncid30);

%% cp90cm
ncid90 = netcdf.create(file_cp90,cmode);

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

vidbio90 = netcdf.defVar(ncid90,'cp90cm','double',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncid90,vidbio90,'CHUNKED',[ni nj 1]);
netcdf.putAtt(ncid90,vidbio90,'long_name','Catch Density of Large Pelagics >=90cm');
netcdf.putAtt(ncid90,vidbio90,'units','g m-2' );
netcdf.defVarFill(ncid90,vidbio90,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid90,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid90,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncid90,varid,'contact','C. Petrik');
netcdf.putAtt(ncid90,varid,'institution','UC San Diego');
netcdf.putAtt(ncid90,varid,'wet weight:C ratio','9:1');
% netcdf.putAtt(ncid90,varid,'feisty fish size','29.24 to 232.08cm');

netcdf.endDef(ncid90);

netcdf.putVar(ncid90,vidlat,lat);
netcdf.putVar(ncid90,vidlon,lon);
netcdf.putVar(ncid90,vidbio90,tpc);
netcdf.putVar(ncid90,vidt90,hist_time);

netcdf.close(ncid90);





