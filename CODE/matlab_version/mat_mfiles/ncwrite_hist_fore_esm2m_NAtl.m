% FEISTY output in N Atl for Andy Visser EU project
% Years 1960-2010 will be good
% annual means for temp & plankton production

clear
close all

%%
gpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/cobalt_data/';
spath = '/Volumes/petrik-lab/Feisty/GCM_Data/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/NAtl/'];

%% times
Hyr = 1861:2005;
Fyr = 2006:2100;

%% annual mean temp & fluxes
load([gpath 'NAtl_ann_means_ESM2M_Hist_Forecast.mat'],'y1','y2',...
    'HTP','HTB','HDet','HNPP','HZM',...
    'FTP','FTB','FDet','FNPP','FZM',...
    'tlat','tlon',...
    'grid','grid_natl','inatl','ID','IDnatl');

%% negs to zero
HZM(HZM(:)<0) = 0.0;
HNPP(HNPP(:)<0) = 0.0;
HDet(HDet(:)<0) = 0.0;

FZM(FZM(:)<0) = 0.0;
FNPP(FNPP(:)<0) = 0.0;
FDet(FDet(:)<0) = 0.0;

%% Quick look
pb = HNPP(:,:,50);
db = HDet(:,:,50);
cb = HZM(:,:,50);

figure(1)
pcolor(log10(pb + eps))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('NPP')

figure(2)
pcolor(log10(db + eps))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('Det')

figure(3)
pcolor(log10(cb + eps))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('Zmeso')

%% netcdf write
% nans to a large number
HZM(isnan(HZM)) = 1.000000020040877e20;
HNPP(isnan(HNPP)) = 1.000000020040877e20;
HDet(isnan(HDet)) = 1.000000020040877e20;
HTB(isnan(HTB)) = 1.000000020040877e20;
HTP(isnan(HTP)) = 1.000000020040877e20;

FZM(isnan(FZM)) = 1.000000020040877e20;
FNPP(isnan(FNPP)) = 1.000000020040877e20;
FDet(isnan(FDet)) = 1.000000020040877e20;
FTB(isnan(FTB)) = 1.000000020040877e20;
FTP(isnan(FTP)) = 1.000000020040877e20;

%%
close all

%% Setup netcdf path to store to
fname1 = ['FEISTY_ESM2M_Historic_1861-2005_'];
fname2 = ['FEISTY_ESM2M_Forecast_2006-2100_'];
fname3 = '.nc';

file_hzm = [fpath fname1 'zmeso' fname3];
file_hnpp = [fpath fname1 'npp' fname3];
file_hdet = [fpath fname1 'det_btm' fname3];
file_htb = [fpath fname1 'temp_btm' fname3];
file_htp = [fpath fname1 'temp_pel' fname3];

file_fzm = [fpath fname2 'zmeso' fname3];
file_fnpp = [fpath fname2 'npp' fname3];
file_fdet = [fpath fname2 'det_btm' fname3];
file_ftb = [fpath fname2 'temp_btm' fname3];
file_ftp = [fpath fname2 'temp_pel' fname3];

[ni,nj,nth] = size(HNPP);
[~,~,ntf] = size(FNPP);

%Use Netcdf4 classic
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

%% Netcdf OUTPUTS =================================================

%% Hzm
ncidZM = netcdf.create(file_hzm,cmode);

time_dim = netcdf.defDim(ncidZM,'time',nth);
lon_dim = netcdf.defDim(ncidZM,'nlon',ni);
lat_dim = netcdf.defDim(ncidZM,'nlat',nj);

vidtZM = netcdf.defVar(ncidZM,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidZM,vidtZM,'long_name','time');
netcdf.putAtt(ncidZM,vidtZM,'standard_name','time');
netcdf.putAtt(ncidZM,vidtZM,'calendar','365_day');
netcdf.putAtt(ncidZM,vidtZM,'axis','T');
netcdf.putAtt(ncidZM,vidtZM,'units','year');

vidlon = netcdf.defVar(ncidZM,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidZM,vidlon,'long_name','longitude');
netcdf.putAtt(ncidZM,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidZM,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidZM,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidZM,vidlat,'long_name','latitude');
netcdf.putAtt(ncidZM,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidZM,vidlat,'axis','Y');

vidbioZM = netcdf.defVar(ncidZM,'zmeso','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidZM,vidbioZM,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name','Mean productivity of mesozooplankton');
netcdf.putAtt(ncidZM,vidbioZM,'units','gC m-2 d-1' );
netcdf.defVarFill(ncidZM,vidbioZM,false,1.000000020040877e20);
netcdf.putAtt(ncidZM,vidbioZM,'missing value',1.000000020040877e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidZM,varid,'_FillValue',1.000000020040877e20);
netcdf.putAtt(ncidZM,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidZM,varid,'institution','UCSD');

netcdf.endDef(ncidZM);

netcdf.putVar(ncidZM,vidlat,tlat);
netcdf.putVar(ncidZM,vidlon,tlon);
netcdf.putVar(ncidZM,vidbioZM,HZM);
netcdf.putVar(ncidZM,vidtZM,Hyr);

netcdf.close(ncidZM);

%%
ncdisp(file_hzm)

%% Hnpp
ncidNPP = netcdf.create(file_hnpp,cmode);

time_dim = netcdf.defDim(ncidNPP,'time',nth);
lon_dim = netcdf.defDim(ncidNPP,'nlon',ni);
lat_dim = netcdf.defDim(ncidNPP,'nlat',nj);

vidtNPP = netcdf.defVar(ncidNPP,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidNPP,vidtNPP,'long_name','time');
netcdf.putAtt(ncidNPP,vidtNPP,'standard_name','time');
netcdf.putAtt(ncidNPP,vidtNPP,'units','year' );
netcdf.putAtt(ncidNPP,vidtNPP,'calendar','365_day');
netcdf.putAtt(ncidNPP,vidtNPP,'axis','T');

vidlon = netcdf.defVar(ncidNPP,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidNPP,vidlon,'long_name','longitude');
netcdf.putAtt(ncidNPP,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidNPP,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidNPP,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidNPP,vidlat,'long_name','latitude');
netcdf.putAtt(ncidNPP,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidNPP,vidlat,'axis','Y');

vidbioNPP = netcdf.defVar(ncidNPP,'npp','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidNPP,vidbioNPP,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidNPP,vidbioNPP,'long_name','Mean net primary productivity');
netcdf.putAtt(ncidNPP,vidbioNPP,'units','gC m-2 d-1' );
netcdf.defVarFill(ncidNPP,vidbioNPP,false,1.000000020040877e20);
netcdf.putAtt(ncidNPP,vidbioNPP,'missing value',1.000000020040877e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidNPP,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidNPP,varid,'_FillValue',1.000000020040877e20);
netcdf.putAtt(ncidNPP,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidNPP,varid,'institution','UCSD');

netcdf.endDef(ncidNPP);

netcdf.putVar(ncidNPP,vidlat,tlat);
netcdf.putVar(ncidNPP,vidlon,tlon);
netcdf.putVar(ncidNPP,vidbioNPP,HNPP);
netcdf.putVar(ncidNPP,vidtNPP,Hyr);

netcdf.close(ncidNPP);

%% Hdb
ncidDet = netcdf.create(file_hdet,cmode);

time_dim = netcdf.defDim(ncidDet,'time',nth);
lon_dim = netcdf.defDim(ncidDet,'nlon',ni);
lat_dim = netcdf.defDim(ncidDet,'nlat',nj);

vidtDet = netcdf.defVar(ncidDet,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidDet,vidtDet,'long_name','time');
netcdf.putAtt(ncidDet,vidtDet,'standard_name','time');
netcdf.putAtt(ncidDet,vidtDet,'calendar','365_day');
netcdf.putAtt(ncidDet,vidtDet,'axis','T');
netcdf.putAtt(ncidDet,vidtDet,'units','year' );

vidlon = netcdf.defVar(ncidDet,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidDet,vidlon,'long_name','longitude');
netcdf.putAtt(ncidDet,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidDet,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidDet,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidDet,vidlat,'long_name','latitude');
netcdf.putAtt(ncidDet,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidDet,vidlat,'axis','Y');

vidbioDet = netcdf.defVar(ncidDet,'det_btm','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidDet,vidbioDet,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidDet,vidbioDet,'long_name','Mean flux of detritus to seafloor');
netcdf.putAtt(ncidDet,vidbioDet,'units','gC m-2 d-1' );
netcdf.defVarFill(ncidDet,vidbioDet,false,1.000000020040877e20);
netcdf.putAtt(ncidDet,vidbioDet,'missing value',1.000000020040877e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidDet,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidDet,varid,'_FillValue',1.000000020040877e20);
netcdf.putAtt(ncidDet,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidDet,varid,'institution','UCSD');

netcdf.endDef(ncidDet);

netcdf.putVar(ncidDet,vidlat,tlat);
netcdf.putVar(ncidDet,vidlon,tlon);
netcdf.putVar(ncidDet,vidbioDet,HDet);
netcdf.putVar(ncidDet,vidtDet,Hyr);

netcdf.close(ncidDet);



%% H TP
ncidTP = netcdf.create(file_htp,cmode);

time_dim = netcdf.defDim(ncidTP,'time',nth);
lon_dim = netcdf.defDim(ncidTP,'nlon',ni);
lat_dim = netcdf.defDim(ncidTP,'nlat',nj);

vidtTP = netcdf.defVar(ncidTP,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidTP,vidtTP,'long_name','time');
netcdf.putAtt(ncidTP,vidtTP,'standard_name','time');
netcdf.putAtt(ncidTP,vidtTP,'units','year' );
netcdf.putAtt(ncidTP,vidtTP,'calendar','365_day');
netcdf.putAtt(ncidTP,vidtTP,'axis','T');

vidlon = netcdf.defVar(ncidTP,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidTP,vidlon,'long_name','longitude');
netcdf.putAtt(ncidTP,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidTP,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidTP,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidTP,vidlat,'long_name','latitude');
netcdf.putAtt(ncidTP,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidTP,vidlat,'axis','Y');

vidbioTP = netcdf.defVar(ncidTP,'temp_pel','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidTP,vidbioTP,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidTP,vidbioTP,'long_name','Mean pelagic (0-100 m) temperature');
netcdf.putAtt(ncidTP,vidbioTP,'units','degC' );
netcdf.defVarFill(ncidTP,vidbioTP,false,1.000000020040877e20);
netcdf.putAtt(ncidTP,vidbioTP,'missing value',1.000000020040877e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidTP,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidTP,varid,'_FillValue',1.000000020040877e20);
netcdf.putAtt(ncidTP,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidTP,varid,'institution','UCSD');

netcdf.endDef(ncidTP);

netcdf.putVar(ncidTP,vidlat,tlat);
netcdf.putVar(ncidTP,vidlon,tlon);
netcdf.putVar(ncidTP,vidbioTP,HTP);
netcdf.putVar(ncidTP,vidtTP,Hyr);

netcdf.close(ncidTP);

%% H TB
ncidTB = netcdf.create(file_htb,cmode);

time_dim = netcdf.defDim(ncidTB,'time',nth);
lon_dim = netcdf.defDim(ncidTB,'nlon',ni);
lat_dim = netcdf.defDim(ncidTB,'nlat',nj);

vidtTB = netcdf.defVar(ncidTB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidTB,vidtTB,'long_name','time');
netcdf.putAtt(ncidTB,vidtTB,'standard_name','time');
netcdf.putAtt(ncidTB,vidtTB,'calendar','365_day');
netcdf.putAtt(ncidTB,vidtTB,'axis','T');
netcdf.putAtt(ncidTB,vidtTB,'units','year' );

vidlon = netcdf.defVar(ncidTB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidTB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidTB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidTB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidTB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidTB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidTB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidTB,vidlat,'axis','Y');

vidbioTB = netcdf.defVar(ncidTB,'temp_btm','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidTB,vidbioTB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidTB,vidbioTB,'long_name','Mean seafloor temperature');
netcdf.putAtt(ncidTB,vidbioTB,'units','degC' );
netcdf.defVarFill(ncidTB,vidbioTB,false,1.000000020040877e20);
netcdf.putAtt(ncidTB,vidbioTB,'missing value',1.000000020040877e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidTB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidTB,varid,'_FillValue',1.000000020040877e20);
netcdf.putAtt(ncidTB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidTB,varid,'institution','UCSD');

netcdf.endDef(ncidTB);

netcdf.putVar(ncidTB,vidlat,tlat);
netcdf.putVar(ncidTB,vidlon,tlon);
netcdf.putVar(ncidTB,vidbioTB,HTB);
netcdf.putVar(ncidTB,vidtTB,Hyr);

netcdf.close(ncidTB);




%% Fzm ------------------------ forecast --------------------------------
ncidZM = netcdf.create(file_fzm,cmode);

time_dim = netcdf.defDim(ncidZM,'time',ntf);
lon_dim = netcdf.defDim(ncidZM,'nlon',ni);
lat_dim = netcdf.defDim(ncidZM,'nlat',nj);

vidtZM = netcdf.defVar(ncidZM,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidZM,vidtZM,'long_name','time');
netcdf.putAtt(ncidZM,vidtZM,'standard_name','time');
netcdf.putAtt(ncidZM,vidtZM,'calendar','365_day');
netcdf.putAtt(ncidZM,vidtZM,'axis','T');
netcdf.putAtt(ncidZM,vidtZM,'units','year');

vidlon = netcdf.defVar(ncidZM,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidZM,vidlon,'long_name','longitude');
netcdf.putAtt(ncidZM,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidZM,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidZM,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidZM,vidlat,'long_name','latitude');
netcdf.putAtt(ncidZM,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidZM,vidlat,'axis','Y');

vidbioZM = netcdf.defVar(ncidZM,'zmeso','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidZM,vidbioZM,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidZM,vidbioZM,'long_name','Mean productivity of mesozooplankton');
netcdf.putAtt(ncidZM,vidbioZM,'units','gC m-2 d-1' );
netcdf.defVarFill(ncidZM,vidbioZM,false,1.000000020040877e20);
netcdf.putAtt(ncidZM,vidbioZM,'missing value',1.000000020040877e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidZM,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidZM,varid,'_FillValue',1.000000020040877e20);
netcdf.putAtt(ncidZM,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidZM,varid,'institution','UCSD');
netcdf.putAtt(ncidZM,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidZM);

%tzm = single(tzm);
netcdf.putVar(ncidZM,vidlat,tlat);
netcdf.putVar(ncidZM,vidlon,tlon);
netcdf.putVar(ncidZM,vidbioZM,FZM);
netcdf.putVar(ncidZM,vidtZM,Fyr);

netcdf.close(ncidZM);

%% Fnpp
ncidNPP = netcdf.create(file_fnpp,cmode);

time_dim = netcdf.defDim(ncidNPP,'time',ntf);
lon_dim = netcdf.defDim(ncidNPP,'nlon',ni);
lat_dim = netcdf.defDim(ncidNPP,'nlat',nj);

vidtNPP = netcdf.defVar(ncidNPP,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidNPP,vidtNPP,'long_name','time');
netcdf.putAtt(ncidNPP,vidtNPP,'standard_name','time');
netcdf.putAtt(ncidNPP,vidtNPP,'units','year' );
netcdf.putAtt(ncidNPP,vidtNPP,'calendar','365_day');
netcdf.putAtt(ncidNPP,vidtNPP,'axis','T');

vidlon = netcdf.defVar(ncidNPP,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidNPP,vidlon,'long_name','longitude');
netcdf.putAtt(ncidNPP,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidNPP,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidNPP,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidNPP,vidlat,'long_name','latitude');
netcdf.putAtt(ncidNPP,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidNPP,vidlat,'axis','Y');

vidbioNPP = netcdf.defVar(ncidNPP,'npp','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidNPP,vidbioNPP,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidNPP,vidbioNPP,'long_name','Mean net primary productivity');
netcdf.putAtt(ncidNPP,vidbioNPP,'units','gC m-2 d-1' );
netcdf.defVarFill(ncidNPP,vidbioNPP,false,1.000000020040877e20);
netcdf.putAtt(ncidNPP,vidbioNPP,'missing value',1.000000020040877e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidNPP,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidNPP,varid,'_FillValue',1.000000020040877e20);
netcdf.putAtt(ncidNPP,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidNPP,varid,'institution','UCSD');
netcdf.putAtt(ncidNPP,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidNPP);

%tnpp = single(tnpp);
netcdf.putVar(ncidNPP,vidlat,tlat);
netcdf.putVar(ncidNPP,vidlon,tlon);
netcdf.putVar(ncidNPP,vidbioNPP,FNPP);
netcdf.putVar(ncidNPP,vidtNPP,Fyr);

netcdf.close(ncidNPP);

%% Fdb
ncidDet = netcdf.create(file_fdet,cmode);

time_dim = netcdf.defDim(ncidDet,'time',ntf);
lon_dim = netcdf.defDim(ncidDet,'nlon',ni);
lat_dim = netcdf.defDim(ncidDet,'nlat',nj);

vidtDet = netcdf.defVar(ncidDet,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidDet,vidtDet,'long_name','time');
netcdf.putAtt(ncidDet,vidtDet,'standard_name','time');
netcdf.putAtt(ncidDet,vidtDet,'calendar','365_day');
netcdf.putAtt(ncidDet,vidtDet,'axis','T');
netcdf.putAtt(ncidDet,vidtDet,'units','year' );

vidlon = netcdf.defVar(ncidDet,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidDet,vidlon,'long_name','longitude');
netcdf.putAtt(ncidDet,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidDet,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidDet,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidDet,vidlat,'long_name','latitude');
netcdf.putAtt(ncidDet,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidDet,vidlat,'axis','Y');

vidbioDet = netcdf.defVar(ncidDet,'det_btm','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidDet,vidbioDet,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidDet,vidbioDet,'long_name','Mean flux of detritus to seafloor');
netcdf.putAtt(ncidDet,vidbioDet,'units','gC m-2 d-1' );
netcdf.defVarFill(ncidDet,vidbioDet,false,1.000000020040877e20);
netcdf.putAtt(ncidDet,vidbioDet,'missing value',1.000000020040877e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidDet,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidDet,varid,'_FillValue',1.000000020040877e20);
netcdf.putAtt(ncidDet,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidDet,varid,'institution','UCSD');
netcdf.putAtt(ncidDet,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidDet);

%tdb = single(tdb);
netcdf.putVar(ncidDet,vidlat,tlat);
netcdf.putVar(ncidDet,vidlon,tlon);
netcdf.putVar(ncidDet,vidbioDet,FDet);
netcdf.putVar(ncidDet,vidtDet,Fyr);

netcdf.close(ncidDet);

%%
ncdisp(file_fdet)

%% F TP
ncidTP = netcdf.create(file_ftp,cmode);

time_dim = netcdf.defDim(ncidTP,'time',ntf);
lon_dim = netcdf.defDim(ncidTP,'nlon',ni);
lat_dim = netcdf.defDim(ncidTP,'nlat',nj);

vidtTP = netcdf.defVar(ncidTP,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidTP,vidtTP,'long_name','time');
netcdf.putAtt(ncidTP,vidtTP,'standard_name','time');
netcdf.putAtt(ncidTP,vidtTP,'units','year' );
netcdf.putAtt(ncidTP,vidtTP,'calendar','365_day');
netcdf.putAtt(ncidTP,vidtTP,'axis','T');

vidlon = netcdf.defVar(ncidTP,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidTP,vidlon,'long_name','longitude');
netcdf.putAtt(ncidTP,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidTP,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidTP,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidTP,vidlat,'long_name','latitude');
netcdf.putAtt(ncidTP,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidTP,vidlat,'axis','Y');

vidbioTP = netcdf.defVar(ncidTP,'temp_pel','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidTP,vidbioTP,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidTP,vidbioTP,'long_name','Mean pelagic (0-100 m) temperature');
netcdf.putAtt(ncidTP,vidbioTP,'units','degC' );
netcdf.defVarFill(ncidTP,vidbioTP,false,1.000000020040877e20);
netcdf.putAtt(ncidTP,vidbioTP,'missing value',1.000000020040877e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidTP,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidTP,varid,'_FillValue',1.000000020040877e20);
netcdf.putAtt(ncidTP,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidTP,varid,'institution','UCSD');

netcdf.endDef(ncidTP);

netcdf.putVar(ncidTP,vidlat,tlat);
netcdf.putVar(ncidTP,vidlon,tlon);
netcdf.putVar(ncidTP,vidbioTP,FTP);
netcdf.putVar(ncidTP,vidtTP,Fyr);

netcdf.close(ncidTP);

%% H TB
ncidTB = netcdf.create(file_ftb,cmode);

time_dim = netcdf.defDim(ncidTB,'time',ntf);
lon_dim = netcdf.defDim(ncidTB,'nlon',ni);
lat_dim = netcdf.defDim(ncidTB,'nlat',nj);

vidtTB = netcdf.defVar(ncidTB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidTB,vidtTB,'long_name','time');
netcdf.putAtt(ncidTB,vidtTB,'standard_name','time');
netcdf.putAtt(ncidTB,vidtTB,'calendar','365_day');
netcdf.putAtt(ncidTB,vidtTB,'axis','T');
netcdf.putAtt(ncidTB,vidtTB,'units','year' );

vidlon = netcdf.defVar(ncidTB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidTB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidTB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidTB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidTB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidTB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidTB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidTB,vidlat,'axis','Y');

vidbioTB = netcdf.defVar(ncidTB,'temp_btm','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidTB,vidbioTB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidTB,vidbioTB,'long_name','Mean seafloor temperature');
netcdf.putAtt(ncidTB,vidbioTB,'units','degC' );
netcdf.defVarFill(ncidTB,vidbioTB,false,1.000000020040877e20);
netcdf.putAtt(ncidTB,vidbioTB,'missing value',1.000000020040877e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidTB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidTB,varid,'_FillValue',1.000000020040877e20);
netcdf.putAtt(ncidTB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidTB,varid,'institution','UCSD');

netcdf.endDef(ncidTB);

netcdf.putVar(ncidTB,vidlat,tlat);
netcdf.putVar(ncidTB,vidlon,tlon);
netcdf.putVar(ncidTB,vidbioTB,FTB);
netcdf.putVar(ncidTB,vidtTB,Fyr);

netcdf.close(ncidTB);



