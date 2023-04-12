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

%% SD
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_v3.2_empHP_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

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

allDPB = SD.bio + MD.bio + LD.bio + allPB;
allDPC = MD.yield + LD.yield + allPC;

clear SD MD LD allPB allPC

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

allFDPB = SF.bio + MF.bio + allDPB;
allTC = MF.yield + allDPC;

clear SF MF allDPB allDPC

%% Benthic material
ncid = netcdf.open([fpath 'Hist_ctrlclim_All_fishobs_v3.2_empHP_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

allTB = biomass + allFDPB;
clear biomass allFDPB

%% catch totals per month
% mult mean catch per day by # days each mo
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
mos = repmat(MNTH,nid,(nt/12));

allTC = allTC .* mos;

%% ========================== netcdf info ===============================

%% ESM
epath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
load([epath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'],'GRD');
load([epath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'],'LAT','LON');

[ni,nj] = size(LAT);

lat = LAT(1,:);
lon = LON(:,1);

%% Reshape to lat,lon,yr
[nid,nt] = size(allTB);

allCB = 1.0e20*ones(ni,nj,nt);
allCC = allCB;

for y=1:nt
    gtcb = 1.0e20*ones(ni,nj);
    ttcb = allTB(:,y);
    gtcb(GRD.ID) = ttcb;
    allCB(:,:,y) = gtcb;
    
    gtcc = 1.0e20*ones(ni,nj);
    ttcc = allTC(:,y);
    gtcc(GRD.ID) = ttcc;
    allCC(:,:,y) = gtcc;
    
end

%%
clear allTB allTC

%% Quick look
pb = allCB(:,:,150);
db = allCC(:,:,150);

figure(1)
pcolor(log10(pb))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('tcb')

figure(2)
pcolor(log10(db))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('tc')

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

file_tcb = [fpath fname1 'tcb' fname2];

file_tcc = [fpath fname1 'tc' fname2];

[ni,nj,nt] = size(allCB);

cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

%% tcb
ncidCB = netcdf.create(file_tcb,cmode);

lon_dim = netcdf.defDim(ncidCB,'lon',ni);
lat_dim = netcdf.defDim(ncidCB,'lat',nj);
time_dim = netcdf.defDim(ncidCB,'time',nt);

vidlat = netcdf.defVar(ncidCB,'lat','double',lat_dim);
netcdf.putAtt(ncidCB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidCB,vidlat,'standard_name','lat');
netcdf.putAtt(ncidCB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidCB,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidCB,'lon','double',lon_dim);
netcdf.putAtt(ncidCB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidCB,vidlon,'standard_name','lon');
netcdf.putAtt(ncidCB,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidCB,vidlon,'axis','X');

vidtCB = netcdf.defVar(ncidCB,'time','double',time_dim);
netcdf.putAtt(ncidCB,vidtCB,'long_name',time_long_name);
netcdf.putAtt(ncidCB,vidtCB,'long_name',time_standard_name);
netcdf.putAtt(ncidCB,vidtCB,'calendar',calendar);
netcdf.putAtt(ncidCB,vidtCB,'axis',time_axis);
netcdf.putAtt(ncidCB,vidtCB,'units',time_units);

vidbioCB = netcdf.defVar(ncidCB,'tcb','double',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidCB,vidbioCB,'CHUNKED',[ni nj 1]);
netcdf.putAtt(ncidCB,vidbioCB,'long_name','Total Consumer Biomass Density');
netcdf.putAtt(ncidCB,vidbioCB,'units','g m-2' );
netcdf.defVarFill(ncidCB,vidbioCB,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidCB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidCB,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidCB,varid,'contact','C. Petrik');
netcdf.putAtt(ncidCB,varid,'institution','UC San Diego');
netcdf.putAtt(ncidCB,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidCB,varid,'includes benthos','yes');

netcdf.endDef(ncidCB);

netcdf.putVar(ncidCB,vidlat,lat);
netcdf.putVar(ncidCB,vidlon,lon);
netcdf.putVar(ncidCB,vidbioCB,allCB);
netcdf.putVar(ncidCB,vidtCB,hist_time);

netcdf.close(ncidCB);


% tc -------------------------CATCH-----------------------------
ncidCB = netcdf.create(file_tcc,cmode);

lon_dim = netcdf.defDim(ncidCB,'lon',ni);
lat_dim = netcdf.defDim(ncidCB,'lat',nj);
time_dim = netcdf.defDim(ncidCB,'time',nt);

vidlat = netcdf.defVar(ncidCB,'lat','double',lat_dim);
netcdf.putAtt(ncidCB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidCB,vidlat,'standard_name','lat');
netcdf.putAtt(ncidCB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidCB,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidCB,'lon','double',lon_dim);
netcdf.putAtt(ncidCB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidCB,vidlon,'standard_name','lon');
netcdf.putAtt(ncidCB,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidCB,vidlon,'axis','X');

vidtCB = netcdf.defVar(ncidCB,'time','double',time_dim);
netcdf.putAtt(ncidCB,vidtCB,'long_name',time_long_name);
netcdf.putAtt(ncidCB,vidtCB,'long_name',time_standard_name);
netcdf.putAtt(ncidCB,vidtCB,'calendar',calendar);
netcdf.putAtt(ncidCB,vidtCB,'axis',time_axis);
netcdf.putAtt(ncidCB,vidtCB,'units',time_units);

vidbioCB = netcdf.defVar(ncidCB,'tc','double',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidCB,vidbioCB,'CHUNKED',[ni nj 1]);
netcdf.putAtt(ncidCB,vidbioCB,'long_name','Total Catch');
netcdf.putAtt(ncidCB,vidbioCB,'units','g m-2' );
netcdf.defVarFill(ncidCB,vidbioCB,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidCB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidCB,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidCB,varid,'contact','C. Petrik');
netcdf.putAtt(ncidCB,varid,'institution','UC San Diego');
netcdf.putAtt(ncidCB,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidCB,varid,'includes benthos','no');

netcdf.endDef(ncidCB);

netcdf.putVar(ncidCB,vidlat,lat);
netcdf.putVar(ncidCB,vidlon,lon);
netcdf.putVar(ncidCB,vidbioCB,allCC);
netcdf.putVar(ncidCB,vidtCB,hist_time);

netcdf.close(ncidCB);


