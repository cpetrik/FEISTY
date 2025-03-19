% FEISTY output at all locations

clear 
close all

%%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

fpath=['/project/Feisty/NC/Matlab_new_size/' cfile '/QuarterDeg/'];

mod = 'Hindcast_All_fish03_v2_';

%% SP
ncid = netcdf.open([fpath mod 'sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%%
[ni,nt] = size(biomass);

SP.bio = biomass;
clear biomass 

%% SF
ncid = netcdf.open([fpath mod 'sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
clear biomass 

% SD
ncid = netcdf.open([fpath mod 'sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
clear biomass

% MP
ncid = netcdf.open([fpath mod 'med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
clear biomass 

% MF
ncid = netcdf.open([fpath mod 'med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
clear biomass 

% MD
ncid = netcdf.open([fpath mod 'med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
clear biomass 

% LP
ncid = netcdf.open([fpath mod 'lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
clear biomass 

% LD
ncid = netcdf.open([fpath mod 'lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
clear biomass 

% Benthic material
ncid = netcdf.open([fpath mod 'bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
clear biomass

%% Netcdf OUTPUTS =================================================
t=time;
mo=t/12;
yr=mo+1948;

allF = SF.bio + MF.bio;
allP = SP.bio + MP.bio + LP.bio;
allD = SD.bio + MD.bio + LD.bio;
allB = Bent.bio;
allS = SF.bio + SP.bio + SD.bio;
allM = MF.bio + MP.bio + MD.bio;
allL = LP.bio + LD.bio;

%% Reshape to lat,lon,yr
AllF = NaN*ones(ni,nj,nt);
AllP = NaN*ones(ni,nj,nt);
AllD = NaN*ones(ni,nj,nt);
AllB = NaN*ones(ni,nj,nt);
AllS = NaN*ones(ni,nj,nt);
AllM = NaN*ones(ni,nj,nt);
AllL = NaN*ones(ni,nj,nt);

for z=1:nt
    Zf=NaN*ones(ni,nj);
    Zp=NaN*ones(ni,nj);
    Zd=NaN*ones(ni,nj);
    Zs=NaN*ones(ni,nj);
    Zm=NaN*ones(ni,nj);
    Zl=NaN*ones(ni,nj);
    Zb=NaN*ones(ni,nj);
    
    Zf(GRD.ID)=allF(:,z);
    Zp(GRD.ID)=allP(:,z);
    Zd(GRD.ID)=allD(:,z);
    Zs(GRD.ID)=allS(:,z);
    Zm(GRD.ID)=allM(:,z);
    Zl(GRD.ID)=allL(:,z);
    Zb(GRD.ID)=allB(:,z);
    
    AllF(:,:,z) = Zf;
    AllP(:,:,z) = Zp;
    AllD(:,:,z) = Zd;
    AllB(:,:,z) = Zb;
    AllS(:,:,z) = Zs;
    AllM(:,:,z) = Zm;
    AllL(:,:,z) = Zl;
end


%%
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
% save([fpath 'CESM_DPLE_CCLME_outputs_monthly_' harv cfile '.mat'],...
%     'time','mo','lat2','lon2','cccol','ccrow',...
%     'AllF','AllD','AllP','AllB',...
%     'AllS','AllM','AllL');


%% netcdf write
% nans to a negative number
AllF(isnan(AllF)) = 1.000000020040877e-20;
AllP(isnan(AllP)) = 1.000000020040877e-20;
AllD(isnan(AllD)) = 1.000000020040877e-20;
AllB(isnan(AllB)) = 1.000000020040877e-20;
AllS(isnan(AllS)) = 1.000000020040877e-20;
AllM(isnan(AllM)) = 1.000000020040877e-20;
AllL(isnan(AllL)) = 1.000000020040877e-20;

AllFish = AllF + AllP + AllD;

%% Quick look
pb = AllP(:,:,50);
db = AllD(:,:,50);
cb = AllF(:,:,50);
bp30 = AllB(:,:,50);

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

%%
close all

%% Setup netcdf path to store to
fname1 = 'feisty_cesm-fosi_monthly_';
fname3 = '.nc';

file_tfb = [fpath fname1 'tfb' fname3];
file_tpb = [fpath fname1 'tpb' fname3];
file_tdb = [fpath fname1 'tdb' fname3];
file_tbb = [fpath fname1 'tbb' fname3];
file_tsb = [fpath fname1 'tsb' fname3];
file_tmb = [fpath fname1 'tmb' fname3];
file_tlb = [fpath fname1 'tlb' fname3];
file_tvb = [fpath fname1 'tvb' fname3];

[ni,nj,nt] = size(AllP);

%%
LAT = TLAT;
LON = TLONG;

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

vidbioFB = netcdf.defVar(ncidFB,'tfb','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidFB,vidbioFB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidFB,vidbioFB,'long_name','Biomass of Forage Fish');
netcdf.putAtt(ncidFB,vidbioFB,'units','gWW m-2' );
netcdf.defVarFill(ncidFB,vidbioFB,false,1.000000020040877e-20);
netcdf.putAtt(ncidFB,vidbioFB,'missing value',1.000000020040877e-20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidFB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidFB,varid,'_FillValue',1.000000020040877e-20);
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

vidbioPB = netcdf.defVar(ncidPB,'tpb','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidPB,vidbioPB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidPB,vidbioPB,'long_name','Biomass of Large Pelagic Fish');
netcdf.putAtt(ncidPB,vidbioPB,'units','gWW m-2' );
netcdf.defVarFill(ncidPB,vidbioPB,false,1.000000020040877e-20);
netcdf.putAtt(ncidPB,vidbioPB,'missing value',1.000000020040877e-20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidPB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidPB,varid,'_FillValue',1.000000020040877e-20);
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

vidbioDB = netcdf.defVar(ncidDB,'tdb','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidDB,vidbioDB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidDB,vidbioDB,'long_name','Biomass of Demersal Fish');
netcdf.putAtt(ncidDB,vidbioDB,'units','gWW m-2' );
netcdf.defVarFill(ncidDB,vidbioDB,false,1.000000020040877e-20);
netcdf.putAtt(ncidDB,vidbioDB,'missing value',1.000000020040877e-20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidDB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidDB,varid,'_FillValue',1.000000020040877e-20);
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

%% tbb
ncidBB = netcdf.create(file_tbb,cmode);

time_dim = netcdf.defDim(ncidBB,'time',nt);
lon_dim = netcdf.defDim(ncidBB,'nlon',ni);
lat_dim = netcdf.defDim(ncidBB,'nlat',nj);

vidtBB = netcdf.defVar(ncidBB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidBB,vidtBB,'long_name','time');
netcdf.putAtt(ncidBB,vidtBB,'standard_name','time');
netcdf.putAtt(ncidBB,vidtBB,'calendar','365_day');
netcdf.putAtt(ncidBB,vidtBB,'axis','T');
netcdf.putAtt(ncidBB,vidtBB,'units','year' );

vidlon = netcdf.defVar(ncidBB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidBB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidBB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidBB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidBB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidBB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidBB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidBB,vidlat,'axis','Y');

vidbioBB = netcdf.defVar(ncidBB,'tbb','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidBB,vidbioBB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidBB,vidbioBB,'long_name','Biomass of Benthic Invertebrates');
netcdf.putAtt(ncidBB,vidbioBB,'units','g  m-2' );
netcdf.defVarFill(ncidBB,vidbioBB,false,1.000000020040877e-20);
netcdf.putAtt(ncidBB,vidbioBB,'missing value',1.000000020040877e-20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidBB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidBB,varid,'_FillValue',1.000000020040877e-20);
netcdf.putAtt(ncidBB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidBB,varid,'institution','UCSD');
netcdf.putAtt(ncidBB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidBB);

%tbb = single(tbb);
netcdf.putVar(ncidBB,vidlat,LAT);
netcdf.putVar(ncidBB,vidlon,LON);
netcdf.putVar(ncidBB,vidbioBB,AllB);
netcdf.putVar(ncidBB,vidtBB,yr);

netcdf.close(ncidBB);

%% tsb
ncidSB = netcdf.create(file_tsb,cmode);

time_dim = netcdf.defDim(ncidSB,'time',nt);
lon_dim = netcdf.defDim(ncidSB,'nlon',ni);
lat_dim = netcdf.defDim(ncidSB,'nlat',nj);

vidtSB = netcdf.defVar(ncidSB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidSB,vidtSB,'long_name','time');
netcdf.putAtt(ncidSB,vidtSB,'standard_name','time');
netcdf.putAtt(ncidSB,vidtSB,'calendar','365_day');
netcdf.putAtt(ncidSB,vidtSB,'axis','T');
netcdf.putAtt(ncidSB,vidtSB,'units','year' );

vidlon = netcdf.defVar(ncidSB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidSB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidSB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidSB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'axis','Y');

vidbioSB = netcdf.defVar(ncidSB,'tsb','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidSB,vidbioSB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidSB,vidbioSB,'long_name','Biomass of Small Fish');
netcdf.putAtt(ncidSB,vidbioSB,'units','gWW m-2' );
netcdf.defVarFill(ncidSB,vidbioSB,false,1.000000020040877e-20);
netcdf.putAtt(ncidSB,vidbioSB,'missing value',1.000000020040877e-20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSB,varid,'_FillValue',1.000000020040877e-20);
netcdf.putAtt(ncidSB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidSB,varid,'institution','UCSD');
netcdf.putAtt(ncidSB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidSB);

%tfb = single(tfb);
netcdf.putVar(ncidSB,vidlat,LAT);
netcdf.putVar(ncidSB,vidlon,LON);
netcdf.putVar(ncidSB,vidbioSB,AllS);
netcdf.putVar(ncidSB,vidtSB,yr);

netcdf.close(ncidSB);

%% tmb
ncidMB = netcdf.create(file_tmb,cmode);

time_dim = netcdf.defDim(ncidMB,'time',nt);
lon_dim = netcdf.defDim(ncidMB,'nlon',ni);
lat_dim = netcdf.defDim(ncidMB,'nlat',nj);

vidtMB = netcdf.defVar(ncidMB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidMB,vidtMB,'long_name','time');
netcdf.putAtt(ncidMB,vidtMB,'standard_name','time');
netcdf.putAtt(ncidMB,vidtMB,'units','year' );
netcdf.putAtt(ncidMB,vidtMB,'calendar','365_day');
netcdf.putAtt(ncidMB,vidtMB,'axis','T');

vidlon = netcdf.defVar(ncidMB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidMB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidMB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidMB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidMB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidMB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidMB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidMB,vidlat,'axis','Y');

vidbioMB = netcdf.defVar(ncidMB,'tmb','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidMB,vidbioMB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidMB,vidbioMB,'long_name','Biomass of Medium Fish');
netcdf.putAtt(ncidMB,vidbioMB,'units','gWW m-2' );
netcdf.defVarFill(ncidMB,vidbioMB,false,1.000000020040877e-20);
netcdf.putAtt(ncidMB,vidbioMB,'missing value',1.000000020040877e-20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidMB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidMB,varid,'_FillValue',1.000000020040877e-20);
netcdf.putAtt(ncidMB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidMB,varid,'institution','UCSD');
netcdf.putAtt(ncidMB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidMB);

%tpb = single(tpb);
netcdf.putVar(ncidMB,vidlat,LAT);
netcdf.putVar(ncidMB,vidlon,LON);
netcdf.putVar(ncidMB,vidbioMB,AllM);
netcdf.putVar(ncidMB,vidtMB,yr);

netcdf.close(ncidMB);

%% tlb
ncidLB = netcdf.create(file_tlb,cmode);

time_dim = netcdf.defDim(ncidLB,'time',nt);
lon_dim = netcdf.defDim(ncidLB,'nlon',ni);
lat_dim = netcdf.defDim(ncidLB,'nlat',nj);

vidtLB = netcdf.defVar(ncidLB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidLB,vidtLB,'long_name','time');
netcdf.putAtt(ncidLB,vidtLB,'standard_name','time');
netcdf.putAtt(ncidLB,vidtLB,'calendar','365_day');
netcdf.putAtt(ncidLB,vidtLB,'axis','T');
netcdf.putAtt(ncidLB,vidtLB,'units','year' );

vidlon = netcdf.defVar(ncidLB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidLB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidLB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidLB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidLB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidLB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidLB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidLB,vidlat,'axis','Y');

vidbioLB = netcdf.defVar(ncidLB,'tlb','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidLB,vidbioLB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidLB,vidbioLB,'long_name','Biomass of Large Fish');
netcdf.putAtt(ncidLB,vidbioLB,'units','gWW m-2' );
netcdf.defVarFill(ncidLB,vidbioLB,false,1.000000020040877e-20);
netcdf.putAtt(ncidLB,vidbioLB,'missing value',1.000000020040877e-20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidLB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidLB,varid,'_FillValue',1.000000020040877e-20);
netcdf.putAtt(ncidLB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidLB,varid,'institution','UCSD');
netcdf.putAtt(ncidLB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidLB);

%tdb = single(tdb);
netcdf.putVar(ncidLB,vidlat,LAT);
netcdf.putVar(ncidLB,vidlon,LON);
netcdf.putVar(ncidLB,vidbioLB,AllL);
netcdf.putVar(ncidLB,vidtLB,yr);

netcdf.close(ncidLB);

%% All fish
ncidVB = netcdf.create(file_tvb,cmode);

time_dim = netcdf.defDim(ncidVB,'time',nt);
lon_dim = netcdf.defDim(ncidVB,'nlon',ni);
lat_dim = netcdf.defDim(ncidVB,'nlat',nj);

vidtVB = netcdf.defVar(ncidVB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidVB,vidtVB,'long_name','time');
netcdf.putAtt(ncidVB,vidtVB,'standard_name','time');
netcdf.putAtt(ncidVB,vidtVB,'calendar','365_day');
netcdf.putAtt(ncidVB,vidtVB,'axis','T');
netcdf.putAtt(ncidVB,vidtVB,'units','year' );

vidlon = netcdf.defVar(ncidVB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidVB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidVB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidVB,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidVB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidVB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidVB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidVB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidVB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidVB,vidlat,'axis','Y');

vidbioVB = netcdf.defVar(ncidVB,'tvb','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidVB,vidbioVB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidVB,vidbioVB,'long_name','Biomass of All Fish');
netcdf.putAtt(ncidVB,vidbioVB,'units','gWW m-2' );
netcdf.defVarFill(ncidVB,vidbioVB,false,1.000000020040877e-20);
netcdf.putAtt(ncidVB,vidbioVB,'missing value',1.000000020040877e-20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidVB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidVB,varid,'_FillValue',1.000000020040877e-20);
netcdf.putAtt(ncidVB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidVB,varid,'institution','UCSD');
netcdf.putAtt(ncidVB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidVB);

%tbb = single(tbb);
netcdf.putVar(ncidVB,vidlat,LAT);
netcdf.putVar(ncidVB,vidlon,LON);
netcdf.putVar(ncidVB,vidbioVB,AllFish);
netcdf.putVar(ncidVB,vidtVB,yr);

netcdf.close(ncidVB);


