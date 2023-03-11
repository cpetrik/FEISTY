% FEISTY output at all locations
% PI ctrlclim pristine 1/4 degree
% Early periods (1841-1960) can have 1841 as the reference year 
% Historic experimental period uses 1901 as the reference year

clear
close all

load('FishMIP_phase3a_exper_times.mat')

time_long_name = 'time';
time_standard_name = 'time';
time_units = 'days since 1841-1-1 00:00:00'; %pi_units; days since xxxx-1-1 00:00:00
time_axis = 'T';
calendar = '365_day';

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];

%% loop over cycles
nid = 670589;
ntc = 240;
nta = ntc*6;
yrs = 1841:1960;
st = 1:ntc:nta;
en = ntc:ntc:nta;

SP.bio = nan*ones(nid,nta);

%%
for c=1:6

    %% SP
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_sml_p_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    %%
    SP.bio(:,st(c):en(c)) = biomass;
    clear biomass

    %% SF
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_sml_f_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    SF.bio(:,st(c):en(c)) = biomass;
    clear biomass

    % SD
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_sml_d_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    SD.bio(:,st(c):en(c)) = biomass;
    clear biomass

    % MP
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_med_p_cycle',num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);...

    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MP.bio(:,st(c):en(c)) = biomass;
    clear biomass

    % MF
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_med_f_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MF.bio(:,st(c):en(c)) = biomass;
    clear biomass

    % MD
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_med_d_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MD.bio(:,st(c):en(c)) = biomass;
    clear biomass

    % LP
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_lrg_p_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    LP.bio(:,st(c):en(c)) = biomass;
    clear biomass

    % LD
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_lrg_d_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    LD.bio(:,st(c):en(c)) = biomass;
    clear biomass

    % Benthic material
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_bent_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    Bent.bio(:,st(c):en(c)) = biomass;
    clear biomass

end

allPB = SP.bio + MP.bio + LP.bio;
allDB = SD.bio + MD.bio + LD.bio;
allFB = SF.bio + MF.bio;
allBB = Bent.bio;

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
tdb = tpb;
tfb = tpb;
tbb = tpb;

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
end

%%
allPelB = tfb + tpb;
allCB = tfb + tpb + tdb + tbb;

%% Quick look
pb = tpb(:,:,150);
db = tdb(:,:,150);
fb = tfb(:,:,150);

figure(1)
pcolor(log10(pb'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('all LPel')

figure(2)
pcolor(log10(db'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('all Dem')

figure(3)
pcolor(log10(fb'))
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
% sens scenario:      15arcmin or 60arcmin

%e.g.
%boats_gfdl-mom6_cobalt2_none_obsclim_histsoc_15arcmin_tcb_global_monthly_1840_2010.nc

close all

%% Setup netcdf path to store to
fname1 = 'feisty_gfdl-mom6-cobalt2_nat_1955-riverine-input_';
fname2 = '_global_monthly_1841_1960.nc';

file_tpb = [fpath fname1 'tpb' fname2];
file_tdb = [fpath fname1 'tdb' fname2];
file_tcb = [fpath fname1 'tcb' fname2];
file_bp30 = [fpath fname1 'bp30cm' fname2];
file_bp90 = [fpath fname1 'bp90cm' fname2];
file_bd90 = [fpath fname1 'bd90cm' fname2];

[ni,nj,nt] = size(tpb);

%% tpb
ncidSB = netcdf.create(file_tpb,'netcdf4');

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
netcdf.putVar(ncidSB,vidtSB,pi_time);

netcdf.close(ncidSB);
%
ncdisp(file_tpb)

%% tdb
ncidSD = netcdf.create(file_tdb,'netcdf4');

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
netcdf.putVar(ncidSD,vidtSD,pi_time);

netcdf.close(ncidSD);

%% tcb
ncidCB = netcdf.create(file_tcb,'netcdf4');

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
netcdf.putVar(ncidCB,vidtCB,pi_time);

netcdf.close(ncidCB);

%% bp30cm
ncid30 = netcdf.create(file_bp30,'netcdf4');

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
netcdf.putVar(ncid30,vidt30,pi_time);

netcdf.close(ncid30);

% bp90cm
ncid90 = netcdf.create(file_bp90,'netcdf4');

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
netcdf.putVar(ncid90,vidt90,pi_time);

netcdf.close(ncid90);

% bd90cm
ncid90 = netcdf.create(file_bd90,'netcdf4');

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
netcdf.putVar(ncid90,vidt90,pi_time);

netcdf.close(ncid90);

