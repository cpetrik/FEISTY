% FEISTY output at all locations
% Spin ctrlclim pristine 1 degree
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

fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];

%% SP
ncid = netcdf.open([fpath 'Spinup_ctrlclim_All_fishobs_v3.2_empHP_sml_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Spinup_ctrlclim_All_fishobs_v3.2_empHP_med_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Spinup_ctrlclim_All_fishobs_v3.2_empHP_lrg_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Spinup_ctrlclim_All_fishobs_v3.2_empHP_sml_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Spinup_ctrlclim_All_fishobs_v3.2_empHP_med_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Spinup_ctrlclim_All_fishobs_v3.2_empHP_lrg_d.nc'],'NC_NOWRITE');
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

%% SF 
ncid = netcdf.open([fpath 'Spinup_ctrlclim_All_fishobs_v3.2_empHP_sml_f.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Spinup_ctrlclim_All_fishobs_v3.2_empHP_med_f.nc'],'NC_NOWRITE');
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

%% Benthic material
ncid = netcdf.open([fpath 'Spinup_ctrlclim_All_fishobs_v3.2_empHP_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

allBB = biomass;
clear biomass

%% catch totals per month
% mult mean catch per day by # days each mo
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
mos = repmat(MNTH,nid,(nt/12));

allFC = allFC .* mos;
allPC = allPC .* mos;
allDC = allDC .* mos;

%% ========================== netcdf info ===============================

%% ESM
epath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';
load([epath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat'],'GRD');
load([epath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat'],'LAT');

[ni,nj] = size(LAT);

lat = LAT(1,:);
lon = LON(:,1);

%% Reshape to lat,lon,yr
[nid,nt] = size(allPB);

tpb = 1.0e20*ones(ni,nj,nt);
tdb = tpb;
tfb = tpb;
tbb = tpb;
tdc = tpb;
tfc = tpb;
tpc = tpb;

for y=1:nt
    gtpb = 1.0e20*ones(ni,nj);
    ttpb = allPB(:,y);
    gtpb(GRD.ID) = ttpb;
    tpb(:,:,y) = gtpb;
    
    gtdb = 1.0e20*ones(ni,nj);
    ttdb = allDB(:,y);
    gtdb(GRD.ID) = ttdb;
    tdb(:,:,y) = gtdb;
    
    gtfb = 1.0e20*ones(ni,nj);
    ttfb = allFB(:,y);
    gtfb(GRD.ID) = ttfb;
    tfb(:,:,y) = gtfb;

    gtbb = 1.0e20*ones(ni,nj);
    ttbb = allBB(:,y);
    gtbb(GRD.ID) = ttbb;
    tbb(:,:,y) = gtbb;

    gtpc = 1.0e20*ones(ni,nj);
    ttpc = allPC(:,y);
    gtpc(GRD.ID) = ttpc;
    tpc(:,:,y) = gtpc;
    
    gtdc = 1.0e20*ones(ni,nj);
    ttdc = allDC(:,y);
    gtdc(GRD.ID) = ttdc;
    tdc(:,:,y) = gtdc;
    
    gtfc = 1.0e20*ones(ni,nj);
    ttfc = allFC(:,y);
    gtfc(GRD.ID) = ttfc;
    tfc(:,:,y) = gtfc;
    
end

%%
allPelB = tfb + tpb;
allCB = tfb + tpb + tdb + tbb;

allPelC = tfc + tpc;
allCC = tfc + tpc + tdc;

allPelB(allPelB>1.0e20) = 1.0e20;
allCB(allCB>1.0e20) = 1.0e20;

allPelC(allPelC>1.0e20) = 1.0e20;
allCC(allCC>1.0e20) = 1.0e20;

%% Quick look
pb = tpb(:,:,150);
db = tdb(:,:,150);
fb = tfb(:,:,150);

figure(1)
pcolor(log10(pb))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('all LPel')

figure(2)
pcolor(log10(db))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('all Dem')

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

% climate scenario:   obsclim or obsclim
% socioecon scenario: histsoc or nat
% sens scenario:      15arcmin or onedeg

%e.g.
%boats_gfdl-mom6_cobalt2_none_obsclim_histsoc_15arcmin_tcb_global_monthly_1840_2010.nc

close all

%% Setup netcdf path to store to
fname1 = 'feisty_gfdl-mom6-cobalt2_obsclim_histsoc_60arcmin_1955-riverine-input_';
fname2 = '_global_monthly_1831_1840.nc';

file_tpb = [fpath fname1 'tpb' fname2];
file_tdb = [fpath fname1 'tdb' fname2];
file_tcb = [fpath fname1 'tcb' fname2];
file_bp30 = [fpath fname1 'bp30cm' fname2];
file_bp90 = [fpath fname1 'bp90cm' fname2];
file_bd90 = [fpath fname1 'bd90cm' fname2];

file_tpc = [fpath fname1 'tpc' fname2];
file_tdc = [fpath fname1 'tdc' fname2];
file_tcc = [fpath fname1 'tc' fname2];
file_cp30 = [fpath fname1 'cp30cm' fname2];
file_cp90 = [fpath fname1 'cp90cm' fname2];
file_cd90 = [fpath fname1 'cd90cm' fname2];

[ni,nj,nt] = size(tpb);

cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

%% tpb
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
netcdf.putVar(ncidSB,vidtSB,spin_time);

netcdf.close(ncidSB);
%%
ncdisp(file_tpb)

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
netcdf.putVar(ncidSD,vidtSD,spin_time);

netcdf.close(ncidSD);

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
netcdf.putVar(ncidCB,vidtCB,spin_time);

netcdf.close(ncidCB);

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
netcdf.putVar(ncid30,vidt30,spin_time);

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
netcdf.putVar(ncid90,vidt90,spin_time);

netcdf.close(ncid90);

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
netcdf.putVar(ncid90,vidt90,spin_time);

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
netcdf.putVar(ncidSB,vidtSB,spin_time);

netcdf.close(ncidSB);

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
netcdf.putVar(ncidSD,vidtSD,spin_time);

netcdf.close(ncidSD);

%% tc
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
netcdf.putVar(ncidCB,vidtCB,spin_time);

netcdf.close(ncidCB);

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
netcdf.putVar(ncid30,vidt30,spin_time);

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
netcdf.putVar(ncid90,vidt90,spin_time);

netcdf.close(ncid90);

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
netcdf.putVar(ncid90,vidt90,spin_time);

netcdf.close(ncid90);



