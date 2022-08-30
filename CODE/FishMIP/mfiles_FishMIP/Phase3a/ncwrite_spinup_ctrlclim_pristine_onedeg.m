% Output of FEISTY after spinup
% 10 cycles of 1961-1980 ctrlclim
% Mean of final month saved to netcdf for regridding fishing effort

clear all
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];

mod = 'pristine_ctrlclim_onedeg';

load([fpath 'Means_Spinup_ctrlclim_pristine_' cfile '.mat'],...
    'mf_mean','mp_mean','md_mean','lp_mean','ld_mean');

% Map data
%cpath = '/Volumes/MIP/Fish-MIP/Phase3/OneDeg/';
cpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';
load([cpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat']);
load([cpath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat'], 'GRD');

[ni,nj]=size(LON);
 
%% Put on grid
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);

Zmf(GRD.ID)=mf_mean;
Zmp(GRD.ID)=mp_mean;
Zmd(GRD.ID)=md_mean;
Zlp(GRD.ID)=lp_mean;
Zld(GRD.ID)=ld_mean;

AllF = Zmf;
AllP = Zmp+Zlp;
AllD = Zmd+Zld;

%% nans to a large number
AllF(isnan(AllF)) = 1.000000020040877e20;
AllP(isnan(AllP)) = 1.000000020040877e20;
AllD(isnan(AllD)) = 1.000000020040877e20;

%% Setup netcdf path to store to
file_tfb = [fpath 'FEISTY_means_spinup_ctrlclim_pristine_onedeg.nc'];

%Use Netcdf4 classic
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

%% tfb
ncidFB = netcdf.create(file_tfb,cmode);

lon_dim = netcdf.defDim(ncidFB,'nlon',ni);
lat_dim = netcdf.defDim(ncidFB,'nlat',nj);

vidlon = netcdf.defVar(ncidFB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidFB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidFB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidFB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidFB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidFB,vidlat,'axis','Y');

vidbioFB = netcdf.defVar(ncidFB,'tfb','NC_FLOAT',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidbioFB,'long_name','Biomass of Forage Fish');
netcdf.putAtt(ncidFB,vidbioFB,'units','gWW m-2' );
netcdf.defVarFill(ncidFB,vidbioFB,false,1.000000020040877e20);
netcdf.putAtt(ncidFB,vidbioFB,'missing value',1.000000020040877e20);

vidbioPB = netcdf.defVar(ncidFB,'tpb','NC_FLOAT',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidbioPB,'long_name','Biomass of Lg Pelagic Fish');
netcdf.putAtt(ncidFB,vidbioPB,'units','gWW m-2' );
netcdf.defVarFill(ncidFB,vidbioPB,false,1.000000020040877e20);
netcdf.putAtt(ncidFB,vidbioPB,'missing value',1.000000020040877e20);

vidbioDB = netcdf.defVar(ncidFB,'tdb','NC_FLOAT',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidbioDB,'long_name','Biomass of Demersal Fish');
netcdf.putAtt(ncidFB,vidbioDB,'units','gWW m-2' );
netcdf.defVarFill(ncidFB,vidbioDB,false,1.000000020040877e20);
netcdf.putAtt(ncidFB,vidbioDB,'missing value',1.000000020040877e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidFB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidFB,varid,'_FillValue',1.000000020040877e20);
netcdf.putAtt(ncidFB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidFB,varid,'institution','UCSD');
netcdf.putAtt(ncidFB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidFB);

netcdf.putVar(ncidFB,vidlat,LAT);
netcdf.putVar(ncidFB,vidlon,LON);
netcdf.putVar(ncidFB,vidbioFB,AllF);
netcdf.putVar(ncidFB,vidbioPB,AllP);
netcdf.putVar(ncidFB,vidbioDB,AllD);

netcdf.close(ncidFB);

%%
ncdisp(file_tfb)


