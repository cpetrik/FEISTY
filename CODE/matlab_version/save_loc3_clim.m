% Save daily forcing of climatology at 3 test locations

clear all
close all

%save to
spath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/';

%! Setup Climatol (loop 5-year climatology of ESM2.6-COBALT)
load('/Volumes/MIP/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_daily.mat');

%! choose where and when to run the model
load('/Volumes/MIP/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_grid.mat');
load('/Users/cpetrik/Dropbox/Princeton/FEISTY_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
locs = [10;13;16];
ID = ids(locs);
names = abbrev(locs);

%%
Tp = COBALT.Tp(ID,:);
Tb = COBALT.Tb(ID,:);
Zm = COBALT.Zm(ID,:) + COBALT.Zl(ID,:);
det= COBALT.det(ID,:);
dZm= COBALT.dZm(ID,:) + COBALT.dZl(ID,:);
H  = GRD.Z(ID);
LAT= GRD.LAT(ID);
LON= GRD.LON(ID);

%% plot
cb=[0/255 50/255 200/255;...    %blue
    238/255 102/255 119/255;... %red
    0 0 0;
    0.50 0.50 0.50];               % grey

set(groot,'defaultAxesColorOrder',cb);

figure
subplot(2,3,1)
plot(1:365,Tp)
title('Tp')

subplot(2,3,2)
plot(1:365,Zm)
title('zooC')

subplot(2,3,3)
plot(1:365,dZm)
title('zoo mort')

subplot(2,3,4)
plot(1:365,Tb)
title('Tb')

subplot(2,3,5)
plot(1:365,det)
title('POC')

subplot(2,3,6)
plot(1:365,det)
title('POC')
ylim([0 0.1])
legend(names)
print('-dpng',[spath 'Clim_locs3_forcing.png'])


%% save to netcdf and mat files
% Setup netcdf path to store to
file_tpb = [spath 'feisty_input_climatol_daily_locs3.nc'];

time = 1:365;
ni = length(ID);
nt = length(time);

%Use Netcdf4 classic
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

%% tpb
ncidSB = netcdf.create(file_tpb,cmode);

time_dim = netcdf.defDim(ncidSB,'time',nt);
loc_dim = netcdf.defDim(ncidSB,'loc',ni);

vidtSB = netcdf.defVar(ncidSB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidSB,vidtSB,'long_name','time');
netcdf.putAtt(ncidSB,vidtSB,'standard_name','time');
netcdf.putAtt(ncidSB,vidtSB,'units','year' );
netcdf.putAtt(ncidSB,vidtSB,'calendar','365_day');
netcdf.putAtt(ncidSB,vidtSB,'axis','T');

vidloc = netcdf.defVar(ncidSB,'loc','NC_DOUBLE',loc_dim);
netcdf.putAtt(ncidSB,vidloc,'long_name','longitude');
netcdf.putAtt(ncidSB,vidloc,'standard_name','longitude');
netcdf.putAtt(ncidSB,vidloc,'units','degrees_east' );
netcdf.putAtt(ncidSB,vidloc,'axis','X');

vidlat = netcdf.defVar(ncidSB,'lat','NC_DOUBLE',loc_dim);
netcdf.putAtt(ncidSB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidSB,vidlat,'axis','Y');

viddep = netcdf.defVar(ncidSB,'dep','NC_DOUBLE',loc_dim);
netcdf.putAtt(ncidSB,viddep,'long_name','bottom depth');
netcdf.putAtt(ncidSB,viddep,'standard_name','bottom depth');
netcdf.putAtt(ncidSB,viddep,'units','m');
netcdf.putAtt(ncidSB,viddep,'axis','Z');

vidTp = netcdf.defVar(ncidSB,'Tp','NC_FLOAT',[loc_dim,time_dim]);
%netcdf.defVarChunking(ncidSB,vidTp,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidSB,vidTp,'long_name','Pelagic mean temperature');
netcdf.putAtt(ncidSB,vidTp,'units','degrees C' );

vidTb = netcdf.defVar(ncidSB,'Tb','NC_FLOAT',[loc_dim,time_dim]);
%netcdf.defVarChunking(ncidSB,vidTb,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidSB,vidTb,'long_name','Bottom temperature');
netcdf.putAtt(ncidSB,vidTb,'units','degrees C' );

vidPOC = netcdf.defVar(ncidSB,'POC','NC_FLOAT',[loc_dim,time_dim]);
%netcdf.defVarChunking(ncidSB,vidPOC,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidSB,vidPOC,'long_name','Particulate organic matter flux to seafloor');
netcdf.putAtt(ncidSB,vidPOC,'units','g m-2 d-1' );

vidzoo = netcdf.defVar(ncidSB,'zoo','NC_FLOAT',[loc_dim,time_dim]);
%netcdf.defVarChunking(ncidSB,vidzoo,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidSB,vidzoo,'long_name','Biomass of mesozooplankton');
netcdf.putAtt(ncidSB,vidzoo,'units','g m-2' );

vidzmort = netcdf.defVar(ncidSB,'zoo_mort','NC_FLOAT',[loc_dim,time_dim]);
%netcdf.defVarChunking(ncidSB,vidzmort,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidSB,vidzmort,'long_name','Mortality loss of mesozooplankton');
netcdf.putAtt(ncidSB,vidzmort,'units','g m-2 d-1' );

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSB,varid,'_FillValue',1.000000020040877e+20);
netcdf.putAtt(ncidSB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidSB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidSB);


netcdf.putVar(ncidSB,vidtSB,time);
netcdf.putVar(ncidSB,vidloc,LON);
netcdf.putVar(ncidSB,vidlat,LAT);
netcdf.putVar(ncidSB,viddep,H);
netcdf.putVar(ncidSB,vidTp,Tp);
netcdf.putVar(ncidSB,vidTb,Tb);
netcdf.putVar(ncidSB,vidPOC,det);
netcdf.putVar(ncidSB,vidzoo,Zm);
netcdf.putVar(ncidSB,vidzmort,dZm);

netcdf.close(ncidSB);
%%
ncdisp(file_tpb)

%%
oGRD = GRD;
GRD.Z = GRD.Z(ID);
GRD.LON = GRD.LON(ID);
GRD.LAT = GRD.LAT(ID);
GRD.lmask = GRD.lmask(ID);
GRD.AREA = GRD.AREA(ID);

%%
ESM.Tp  = COBALT.Tp(ID,:);
ESM.Tb  = COBALT.Tb(ID,:);
ESM.Zm  = COBALT.Zm(ID,:)+COBALT.Zl(ID,:);
ESM.det = COBALT.det(ID,:);
ESM.dZm = COBALT.dZm(ID,:)+COBALT.dZl(ID,:);

%%
save([spath 'feisty_input_climatol_daily_locs3.mat'],'GRD','ESM');

