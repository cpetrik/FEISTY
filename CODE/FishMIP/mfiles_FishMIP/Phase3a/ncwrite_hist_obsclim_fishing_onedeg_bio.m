% FEISTY output at all locations

clear all
close all

%/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/.../OneDeg
%Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100

%/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/.../OneDeg
%Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];

% MP
ncid = netcdf.open([fpath 'Hist_obsclim_All_fishobs_empHP_med_p.nc'],'NC_NOWRITE');
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

% MF
ncid = netcdf.open([fpath 'Hist_obsclim_All_fishobs_empHP_med_f.nc'],'NC_NOWRITE');
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

% MD
ncid = netcdf.open([fpath 'Hist_obsclim_All_fishobs_empHP_med_d.nc'],'NC_NOWRITE');
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

% LP
ncid = netcdf.open([fpath 'Hist_obsclim_All_fishobs_empHP_lrg_p.nc'],'NC_NOWRITE');
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

% LD
ncid = netcdf.open([fpath 'Hist_obsclim_All_fishobs_empHP_lrg_d.nc'],'NC_NOWRITE');
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

%% Benthic material
ncid = netcdf.open([fpath 'Hist_obsclim_All_fishobs_empHP_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
clear biomass

%% Each year
st=1:12:length(time);
en=12:12:length(time);

for n=1:length(st)
    mp_my(:,n)=nanmean(MP.yield(:,st(n):en(n)),2);
    mf_my(:,n)=nanmean(MF.yield(:,st(n):en(n)),2);
    md_my(:,n)=nanmean(MD.yield(:,st(n):en(n)),2);
    lp_my(:,n)=nanmean(LP.yield(:,st(n):en(n)),2);
    ld_my(:,n)=nanmean(LD.yield(:,st(n):en(n)),2);
end

%% Netcdf OUTPUTS =================================================

% Map data
%cpath = '/Volumes/MIP/Fish-MIP/Phase3/OneDeg/';
cpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';
load([cpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat']);
load([cpath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat'], 'GRD');

[ni,nj]=size(LON);
yr=1961:2010;
nt = length(yr);

allF = mf_my;
allP = mp_my + lp_my;
allD = md_my + ld_my;

%% Reshape to lat,lon,yr
AllF = NaN*ones(ni,nj,nt);
AllP = NaN*ones(ni,nj,nt);
AllD = NaN*ones(ni,nj,nt);

for z=1:nt
    Zf=NaN*ones(ni,nj);
    Zp=NaN*ones(ni,nj);
    Zd=NaN*ones(ni,nj);

    Zf(GRD.ID)=allF(:,z);
    Zp(GRD.ID)=allP(:,z);
    Zd(GRD.ID)=allD(:,z);

    AllF(:,:,z) = Zf;
    AllP(:,:,z) = Zp;
    AllD(:,:,z) = Zd;
end


%% netcdf write
% nans to a large number
AllF(isnan(AllF)) = 1.000000020040877e20;
AllP(isnan(AllP)) = 1.000000020040877e20;
AllD(isnan(AllD)) = 1.000000020040877e20;


%% Quick look
pb = AllP(:,:,50);
db = AllD(:,:,50);
cb = AllF(:,:,50);

figure(1)
pcolor(log10(pb'))
shading flat
colormap('jet')
colorbar
caxis([-5 0])
title('allP')

figure(2)
pcolor(log10(db'))
shading flat
colormap('jet')
colorbar
caxis([-5 0])
title('allD')

figure(3)
pcolor(log10(cb'))
shading flat
colormap('jet')
colorbar
caxis([-5 0])
title('all F')

%%
close all

%% Setup netcdf path to store to
fname1 = 'FEISTY_obsclim_fishobs_onedeg_1961-2010_';
fname3 = '.nc';

file_tfb = [fpath fname1 'tfc' fname3];
file_tpb = [fpath fname1 'tpc' fname3];
file_tdb = [fpath fname1 'tdc' fname3];

[ni,nj,nt] = size(AllP);

%Use Netcdf4 classic
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

%% tfb
ncidFB = netcdf.create(file_tfb,cmode);

time_dim = netcdf.defDim(ncidFB,'time',nt);
lon_dim = netcdf.defDim(ncidFB,'nlon',ni);
lat_dim = netcdf.defDim(ncidFB,'nlat',nj);

vidtFB = netcdf.defVar(ncidFB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidFB,vidtFB,'long_name','time');
netcdf.putAtt(ncidFB,vidtFB,'standard_name','time');
netcdf.putAtt(ncidFB,vidtFB,'calendar','365_day');
netcdf.putAtt(ncidFB,vidtFB,'axis','T');
netcdf.putAtt(ncidFB,vidtFB,'units','year' );

vidlon = netcdf.defVar(ncidFB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidFB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidFB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidFB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidFB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidFB,vidlat,'axis','Y');

vidbioFB = netcdf.defVar(ncidFB,'tfc','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidFB,vidbioFB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidFB,vidbioFB,'long_name','Catch of Forage Fish');
netcdf.putAtt(ncidFB,vidbioFB,'units','gWW m-2 d-1' );
netcdf.defVarFill(ncidFB,vidbioFB,false,1.000000020040877e20);
netcdf.putAtt(ncidFB,vidbioFB,'missing value',1.000000020040877e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidFB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidFB,varid,'_FillValue',1.000000020040877e20);
netcdf.putAtt(ncidFB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidFB,varid,'institution','UCSD');
netcdf.putAtt(ncidFB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidFB);

%tfb = single(tfb);
netcdf.putVar(ncidFB,vidlat,LAT);
netcdf.putVar(ncidFB,vidlon,LON);
netcdf.putVar(ncidFB,vidbioFB,AllF);
netcdf.putVar(ncidFB,vidtFB,yr);

netcdf.close(ncidFB);

%%
ncdisp(file_tfb)

%% tpb
ncidPB = netcdf.create(file_tpb,cmode);

time_dim = netcdf.defDim(ncidPB,'time',nt);
lon_dim = netcdf.defDim(ncidPB,'nlon',ni);
lat_dim = netcdf.defDim(ncidPB,'nlat',nj);

vidtPB = netcdf.defVar(ncidPB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidPB,vidtPB,'long_name','time');
netcdf.putAtt(ncidPB,vidtPB,'standard_name','time');
netcdf.putAtt(ncidPB,vidtPB,'units','year' );
netcdf.putAtt(ncidPB,vidtPB,'calendar','365_day');
netcdf.putAtt(ncidPB,vidtPB,'axis','T');

vidlon = netcdf.defVar(ncidPB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidPB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidPB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidPB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidPB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidPB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidPB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidPB,vidlat,'axis','Y');

vidbioPB = netcdf.defVar(ncidPB,'tpc','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidPB,vidbioPB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidPB,vidbioPB,'long_name','Catch of Large Pelagic Fish');
netcdf.putAtt(ncidPB,vidbioPB,'units','gWW m-2 d-1' );
netcdf.defVarFill(ncidPB,vidbioPB,false,1.000000020040877e20);
netcdf.putAtt(ncidPB,vidbioPB,'missing value',1.000000020040877e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidPB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidPB,varid,'_FillValue',1.000000020040877e20);
netcdf.putAtt(ncidPB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidPB,varid,'institution','UCSD');
netcdf.putAtt(ncidPB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidPB);

%tpb = single(tpb);
netcdf.putVar(ncidPB,vidlat,LAT);
netcdf.putVar(ncidPB,vidlon,LON);
netcdf.putVar(ncidPB,vidbioPB,AllP);
netcdf.putVar(ncidPB,vidtPB,yr);

netcdf.close(ncidPB);

%% tdb
ncidDB = netcdf.create(file_tdb,cmode);

time_dim = netcdf.defDim(ncidDB,'time',nt);
lon_dim = netcdf.defDim(ncidDB,'nlon',ni);
lat_dim = netcdf.defDim(ncidDB,'nlat',nj);

vidtDB = netcdf.defVar(ncidDB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidDB,vidtDB,'long_name','time');
netcdf.putAtt(ncidDB,vidtDB,'standard_name','time');
netcdf.putAtt(ncidDB,vidtDB,'calendar','365_day');
netcdf.putAtt(ncidDB,vidtDB,'axis','T');
netcdf.putAtt(ncidDB,vidtDB,'units','year' );

vidlon = netcdf.defVar(ncidDB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidDB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidDB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidDB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidDB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidDB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidDB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidDB,vidlat,'axis','Y');

vidbioDB = netcdf.defVar(ncidDB,'tdc','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidDB,vidbioDB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidDB,vidbioDB,'long_name','Catch of Demersal Fish');
netcdf.putAtt(ncidDB,vidbioDB,'units','gWW m-2 d-1' );
netcdf.defVarFill(ncidDB,vidbioDB,false,1.000000020040877e20);
netcdf.putAtt(ncidDB,vidbioDB,'missing value',1.000000020040877e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidDB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidDB,varid,'_FillValue',1.000000020040877e20);
netcdf.putAtt(ncidDB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidDB,varid,'institution','UCSD');
netcdf.putAtt(ncidDB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidDB);

%tdb = single(tdb);
netcdf.putVar(ncidDB,vidlat,LAT);
netcdf.putVar(ncidDB,vidlon,LON);
netcdf.putVar(ncidDB,vidbioDB,AllD);
netcdf.putVar(ncidDB,vidtDB,yr);

netcdf.close(ncidDB);
