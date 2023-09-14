% FEISTY output in N Atl for Andy Visser EU project
% Years 1960-2010 will be good
% annual means for prod 

clear
close all

%%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/NAtl/'];

%% times
Hyr = 1861:2005;
Fyr = 2006:2100;

%% annual mean biomasses
load([fpath 'NAtl_Means_prod_Historic_' harv '_' cfile '.mat'],...
    'Fall','Pall','Dall');

HFB = Fall;
HPB = Pall;
HDB = Dall;

clear Fall Pall Dall

load([fpath 'NAtl_Means_prod_Forecast_' harv '_' cfile '.mat'],...
    'Fall','Pall','Dall','tlat','tlon');

FFB = Fall;
FPB = Pall;
FDB = Dall;

clear Fall Pall Dall 

%% negs to zero
HFB(HFB(:)<0) = 0.0;
HPB(HPB(:)<0) = 0.0;
HDB(HDB(:)<0) = 0.0;

FFB(FFB(:)<0) = 0.0;
FPB(FPB(:)<0) = 0.0;
FDB(FDB(:)<0) = 0.0;

%% Quick look
pb = HPB(:,:,50);
db = HDB(:,:,50);
cb = HFB(:,:,50);

figure(1)
pcolor(log10(pb + eps))
shading flat
colormap('jet')
colorbar
caxis([-4 0])
title('allP')

figure(2)
pcolor(log10(db + eps))
shading flat
colormap('jet')
colorbar
caxis([-4 0])
title('allD')

figure(3)
pcolor(log10(cb + eps))
shading flat
colormap('jet')
colorbar
caxis([-4 0])
title('all F')

%% netcdf write
% nans to a large number
HFB(isnan(HFB)) = 1.000000020040877e20;
HPB(isnan(HPB)) = 1.000000020040877e20;
HDB(isnan(HDB)) = 1.000000020040877e20;

FFB(isnan(FFB)) = 1.000000020040877e20;
FPB(isnan(FPB)) = 1.000000020040877e20;
FDB(isnan(FDB)) = 1.000000020040877e20;

%%
close all

%% Setup netcdf path to store to
fname1 = ['FEISTY_ESM2M_Historic_1861-2005_'];
fname2 = ['FEISTY_ESM2M_Forecast_2006-2100_'];
fname3 = '.nc';

file_hfb = [fpath fname1 'mfp' fname3];
file_hpb = [fpath fname1 'mpp' fname3];
file_hdb = [fpath fname1 'mdp' fname3];

file_ffb = [fpath fname2 'mfp' fname3];
file_fpb = [fpath fname2 'mpp' fname3];
file_fdb = [fpath fname2 'mdp' fname3];

[ni,nj,nth] = size(HPB);
[~,~,ntf] = size(FPB);

%Use Netcdf4 classic
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

%% Netcdf OUTPUTS =================================================

%% Hfb
ncidFB = netcdf.create(file_hfb,cmode);

time_dim = netcdf.defDim(ncidFB,'time',nth);
lon_dim = netcdf.defDim(ncidFB,'nlon',ni);
lat_dim = netcdf.defDim(ncidFB,'nlat',nj);

vidtFB = netcdf.defVar(ncidFB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidFB,vidtFB,'long_name','time');
netcdf.putAtt(ncidFB,vidtFB,'standard_name','time');
netcdf.putAtt(ncidFB,vidtFB,'calendar','365_day');
netcdf.putAtt(ncidFB,vidtFB,'axis','T');
netcdf.putAtt(ncidFB,vidtFB,'units','year');

vidlon = netcdf.defVar(ncidFB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidFB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidFB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidFB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidFB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidFB,vidlat,'axis','Y');

vidbioFB = netcdf.defVar(ncidFB,'mfp','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidFB,vidbioFB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidFB,vidbioFB,'long_name','Mean productivity of Forage Fish');
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
netcdf.putVar(ncidFB,vidlat,tlat);
netcdf.putVar(ncidFB,vidlon,tlon);
netcdf.putVar(ncidFB,vidbioFB,HFB);
netcdf.putVar(ncidFB,vidtFB,Hyr);

netcdf.close(ncidFB);

%%
ncdisp(file_hfb)

%% Hpb
ncidPB = netcdf.create(file_hpb,cmode);

time_dim = netcdf.defDim(ncidPB,'time',nth);
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

vidbioPB = netcdf.defVar(ncidPB,'mpp','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidPB,vidbioPB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidPB,vidbioPB,'long_name','Mean productivity of Large Pelagic Fish');
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
netcdf.putVar(ncidPB,vidlat,tlat);
netcdf.putVar(ncidPB,vidlon,tlon);
netcdf.putVar(ncidPB,vidbioPB,HPB);
netcdf.putVar(ncidPB,vidtPB,Hyr);

netcdf.close(ncidPB);

%% Hdb
ncidDB = netcdf.create(file_hdb,cmode);

time_dim = netcdf.defDim(ncidDB,'time',nth);
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

vidbioDB = netcdf.defVar(ncidDB,'mdp','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidDB,vidbioDB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidDB,vidbioDB,'long_name','Mean productivity of Demersal Fish');
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
netcdf.putVar(ncidDB,vidlat,tlat);
netcdf.putVar(ncidDB,vidlon,tlon);
netcdf.putVar(ncidDB,vidbioDB,HDB);
netcdf.putVar(ncidDB,vidtDB,Hyr);

netcdf.close(ncidDB);


%% Ffb ------------------------ forecast --------------------------------
ncidFB = netcdf.create(file_ffb,cmode);

time_dim = netcdf.defDim(ncidFB,'time',ntf);
lon_dim = netcdf.defDim(ncidFB,'nlon',ni);
lat_dim = netcdf.defDim(ncidFB,'nlat',nj);

vidtFB = netcdf.defVar(ncidFB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidFB,vidtFB,'long_name','time');
netcdf.putAtt(ncidFB,vidtFB,'standard_name','time');
netcdf.putAtt(ncidFB,vidtFB,'calendar','365_day');
netcdf.putAtt(ncidFB,vidtFB,'axis','T');
netcdf.putAtt(ncidFB,vidtFB,'units','year');

vidlon = netcdf.defVar(ncidFB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidFB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidFB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidFB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidFB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidFB,vidlat,'axis','Y');

vidbioFB = netcdf.defVar(ncidFB,'mfp','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidFB,vidbioFB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidFB,vidbioFB,'long_name','Mean productivity of Forage Fish');
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
netcdf.putVar(ncidFB,vidlat,tlat);
netcdf.putVar(ncidFB,vidlon,tlon);
netcdf.putVar(ncidFB,vidbioFB,FFB);
netcdf.putVar(ncidFB,vidtFB,Fyr);

netcdf.close(ncidFB);

%% Fpb
ncidPB = netcdf.create(file_fpb,cmode);

time_dim = netcdf.defDim(ncidPB,'time',ntf);
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

vidbioPB = netcdf.defVar(ncidPB,'mpp','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidPB,vidbioPB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidPB,vidbioPB,'long_name','Mean productivity of Large Pelagic Fish');
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
netcdf.putVar(ncidPB,vidlat,tlat);
netcdf.putVar(ncidPB,vidlon,tlon);
netcdf.putVar(ncidPB,vidbioPB,FPB);
netcdf.putVar(ncidPB,vidtPB,Fyr);

netcdf.close(ncidPB);

%% Fdb
ncidDB = netcdf.create(file_fdb,cmode);

time_dim = netcdf.defDim(ncidDB,'time',ntf);
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

vidbioDB = netcdf.defVar(ncidDB,'mdp','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidDB,vidbioDB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidDB,vidbioDB,'long_name','Mean productivity of Demersal Fish');
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
netcdf.putVar(ncidDB,vidlat,tlat);
netcdf.putVar(ncidDB,vidlon,tlon);
netcdf.putVar(ncidDB,vidbioDB,FDB);
netcdf.putVar(ncidDB,vidtDB,Fyr);

netcdf.close(ncidDB);

%%
ncdisp(file_fdb)


