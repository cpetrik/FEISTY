% Save time series output of FEISTY forced by CMIP6
% tcb and tc production for Vianney's paper

clear 
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP/FishMIP6/';
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%% gfdl hist
gpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_CMIP6/' cfile '/'];
load([gpath 'Hist_empHP_fishMIP_outputs_monthly_' cfile '.mat'],'time','mo',...
    'allC','allCprod');

% ANNUAL MEAN
[ni,nt] = size(allC);
nyr = length(time)/12;
st=1:12:length(time);
en=12:12:length(time);

GHtcb = nan*ones(ni,nyr);

for n=1:length(st)
    % mean prod
    GHtcb(:,n)=nanmean(allC(:,st(n):en(n)),2);
    
end

GHtcp = allCprod;

clear allC allCprod time mo

%% gfdl ssp585
load([gpath 'SSP585_empHP_fishMIP_outputs_monthly_' cfile '.mat'],'time','mo',...
    'allC','allCprod');

% ANNUAL MEAN
[ni,nt] = size(allC);
nyr = length(time)/12;
st=1:12:length(time);
en=12:12:length(time);

GStcb = nan*ones(ni,nyr);

for n=1:length(st)
    % mean prod
    GStcb(:,n)=nanmean(allC(:,st(n):en(n)),2);
    
end

GStcp = allCprod;

clear allC allCprod time mo

%% ipsl hist
ipath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/IPSL_CMIP6/' cfile '/'];
load([ipath 'Hist_empHP_fishMIP_outputs_monthly_' cfile '.mat'],'time','mo',...
    'allC','allCprod');

% ANNUAL MEAN
[ni,nt] = size(allC);
nyr = length(time)/12;
st=1:12:length(time);
en=12:12:length(time);

IHtcb = nan*ones(ni,nyr);

for n=1:length(st)
    % mean prod
    IHtcb(:,n)=nanmean(allC(:,st(n):en(n)),2);
    
end

IHtcp = allCprod;

clear allC allCprod time mo

%%
load([ipath 'SSP585_empHP_fishMIP_outputs_monthly_' cfile '.mat'],'time','mo',...
    'allC','allCprod');

% ANNUAL MEAN
[ni,nt] = size(allC);
nyr = length(time)/12;
st=1:12:length(time);
en=12:12:length(time);

IStcb = nan*ones(ni,nyr);

for n=1:length(st)
    % mean prod
    IStcb(:,n)=nanmean(allC(:,st(n):en(n)),2);
    
end

IStcp = allCprod;

clear allC allCprod time mo

%% time
yH = 1951:2015;
yS = 2016:2101;

Gtcp = [GHtcp GStcp];
Gtcb = [GHtcb GStcb];
Itcp = [IHtcp IStcp];
Itcb = [IHtcb IStcb];

yr = [yH yS];

%% percent diff

mGbm = mean(Gtcb,'omitnan');
mGpm = mean(Gtcp,'omitnan');
mIbm = mean(Itcb,'omitnan');
mIpm = mean(Itcp,'omitnan');

tid = find(yr>1995 & yr<=2014);

GHbm = mean(mGbm(:,tid));
GHpm = mean(mGpm(:,tid));
IHbm = mean(mIbm(:,tid));
IHpm = mean(mIpm(:,tid));

pdGHtcb = (mGbm-GHbm)./GHbm;
pdGHtcp = (mGpm-GHpm)./GHpm;
pdIHtcb = (mIbm-IHbm)./IHbm;
pdIHtcp = (mIpm-IHpm)./IHpm;

%% viz percent diff from 1995-2014
figure(1)
subplot(2,2,1)
plot(yr,100*pdGHtcp,'r','LineWidth',1); hold on;
plot(yr,100*pdGHtcb,'b','LineWidth',1); hold on;
title('All fish consumers GFDL-FEISTY')
ylabel('% difference from 1995-2014')
legend('production','biomass')
xlim([1950 2100])

subplot(2,2,2)
plot(yr,100*pdIHtcp,'r','LineWidth',1); hold on;
plot(yr,100*pdIHtcb,'b','LineWidth',1); hold on;
title('All fish consumers IPSL-FEISTY')
ylabel('% difference from 1995-2014')
legend('production','biomass')
xlim([1950 2100])

stamp('')
%print('-dpng',[ppath 'Hist_SSP585_ts_empHP_tcb_prod_pdiff_1994_2014.png'])

%% NETCDF SAVE ======================================================

%% Reshape GFDL to lat,lon,yr

load('/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/GFDL/Data_grid_gfdl.mat','GRD');
load('/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/GFDL/gridspec_gfdl_cmip6.mat','LAT','LON');

GID = GRD.ID;

[gni,gnj] = size(LAT);

% Lat & Lon should be vectors
GLAT = LAT(1,:);
GLON = LON(:,1);

clear GRD LAT LON

[nid,nt] = size(Gtcb);

tcbG = 1.0e20*ones(gni,gnj,nt);
tcpG = tcbG;

for y=1:nt
    gtcb = 1.0e20*ones(gni,gnj);
    ttcb = Gtcb(:,y);
    gtcb(GID) = ttcb;
    tcbG(:,:,y) = gtcb;
    
    gtcp = 1.0e20*ones(gni,gnj);
    ttcp = Gtcp(:,y);
    gtcp(GID) = ttcp;
    tcpG(:,:,y) = gtcp;
end

%% Reshape IPSL to lat,lon,yr

load('/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/Data_grid_ipsl.mat','GRD');
load('/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/gridspec_ipsl_cmip6.mat','LAT','LON');

IID = GRD.ID;

[ini,inj] = size(LAT);

% Lat & Lon should be vectors
ILAT = LAT(1,:);
ILON = LON(:,1);

clear GRD LAT LON

[nid,nt] = size(Itcb);

tcbI = 1.0e20*ones(ini,inj,nt);
tcpI = tcbI;

for y=1:nt
    itcb = 1.0e20*ones(ini,inj);
    ttcb = Itcb(:,y);
    itcb(IID) = ttcb;
    tcbI(:,:,y) = itcb;
    
    itcp = 1.0e20*ones(ini,inj);
    ttcp = Itcp(:,y);
    itcp(IID) = ttcp;
    tcpI(:,:,y) = itcp;
end

%% Quick look
Gcp = tcpG(:,:,15);
Gcb = tcbG(:,:,15);
Icp = tcpI(:,:,15);
Icb = tcbI(:,:,15);

figure(2)
pcolor((Gcp'))
shading flat
colormap('jet')
colorbar
caxis([0 100])
title('G prod')

figure(3)
pcolor(log10(Gcb'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('G biom')

figure(4)
pcolor((Icp'))
shading flat
colormap('jet')
colorbar
caxis([0 100])
title('I prod')

figure(5)
pcolor(log10(Icb'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])
title('I biom')

%% Setup netcdf path to store to
close all

fname1 = 'feisty_gfdl_historical_ssp585_nat_default_';
fname2 = '_global_annual_1950_2100.nc';

file_tcbG = [gpath fname1 'tcb' fname2];
file_tcpG = [gpath fname1 'tcp' fname2];

fname3 = 'feisty_ipsl_historical_ssp585_nat_default_';
fname4 = '_global_annual_1950_2100.nc';

file_tcbI = [ipath fname3 'tcb' fname4];
file_tcpI = [ipath fname3 'tcp' fname4];

%% tcb GFDL
ncidCB = netcdf.create(file_tcbG,'netcdf4');

lon_dim = netcdf.defDim(ncidCB,'lon',gni);
lat_dim = netcdf.defDim(ncidCB,'lat',gnj);
time_dim = netcdf.defDim(ncidCB,'time',nt);

vidlat = netcdf.defVar(ncidCB,'lat','double',lat_dim);
netcdf.putAtt(ncidCB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidCB,vidlat,'standard_name','lat');
netcdf.putAtt(ncidCB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidCB,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidCB,'lon','double',lon_dim);
netcdf.putAtt(ncidCB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidCB,vidlon,'standard_name','lon');
netcdf.putAtt(ncidCB,vidlon,'units','degrees_east');
netcdf.putAtt(ncidCB,vidlon,'axis','X');

vidtCB = netcdf.defVar(ncidCB,'time','double',time_dim);
netcdf.putAtt(ncidCB,vidtCB,'long_name','time');
netcdf.putAtt(ncidCB,vidtCB,'standard_name','time');
netcdf.putAtt(ncidCB,vidtCB,'calendar','365_day');
netcdf.putAtt(ncidCB,vidtCB,'axis','T');
netcdf.putAtt(ncidCB,vidtCB,'units','year');

vidbioCB = netcdf.defVar(ncidCB,'tcb','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidCB,vidbioCB,'long_name','total consumer biomass');
netcdf.putAtt(ncidCB,vidbioCB,'units','grams wet weight m-2' );
netcdf.defVarFill(ncidCB,vidbioCB,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidCB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidCB,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidCB,varid,'contact','C. Petrik');
netcdf.putAtt(ncidCB,varid,'institution','UC San Diego');
netcdf.putAtt(ncidCB,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidCB,varid,'includes benthos','yes');

netcdf.endDef(ncidCB);

netcdf.putVar(ncidCB,vidlat,GLAT);
netcdf.putVar(ncidCB,vidlon,GLON);
netcdf.putVar(ncidCB,vidbioCB,tcbG);
netcdf.putVar(ncidCB,vidtCB,yr);

netcdf.close(ncidCB);

%% tcp GFDL
ncidSB = netcdf.create(file_tcpG,'netcdf4');

lon_dim = netcdf.defDim(ncidSB,'lon',gni);
lat_dim = netcdf.defDim(ncidSB,'lat',gnj);
time_dim = netcdf.defDim(ncidSB,'time',nt);

vidlat = netcdf.defVar(ncidSB,'lat','double',lat_dim);
netcdf.putAtt(ncidSB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'standard_name','lat');
netcdf.putAtt(ncidSB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidSB,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidSB,'lon','double',lon_dim);
netcdf.putAtt(ncidSB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidSB,vidlon,'standard_name','lon');
netcdf.putAtt(ncidSB,vidlon,'units','degrees_east');
netcdf.putAtt(ncidSB,vidlon,'axis','X');

vidtSB = netcdf.defVar(ncidSB,'time','double',time_dim);
netcdf.putAtt(ncidSB,vidtSB,'long_name','time');
netcdf.putAtt(ncidSB,vidtSB,'standard_name','time');
netcdf.putAtt(ncidSB,vidtSB,'units','year');
netcdf.putAtt(ncidSB,vidtSB,'calendar','365_day');
netcdf.putAtt(ncidSB,vidtSB,'axis','T');

vidbioSB = netcdf.defVar(ncidSB,'tpb','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSB,vidbioSB,'long_name','total consumer productivity');
netcdf.putAtt(ncidSB,vidbioSB,'units','grams wet weight m-2 y-1' );
netcdf.defVarFill(ncidSB,vidbioSB,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSB,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidSB,varid,'contact','C. Petrik');
netcdf.putAtt(ncidSB,varid,'institution','UC San Diego');
netcdf.putAtt(ncidSB,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidSB,varid,'includes benthos','no');

netcdf.endDef(ncidSB);

netcdf.putVar(ncidSB,vidlat,GLAT);
netcdf.putVar(ncidSB,vidlon,GLON);
netcdf.putVar(ncidSB,vidbioSB,tcpG);
netcdf.putVar(ncidSB,vidtSB,yr);

netcdf.close(ncidSB);

%%
ncdisp(file_tcpG)

%% tcb IPSL
ncidCB = netcdf.create(file_tcbI,'netcdf4');

lon_dim = netcdf.defDim(ncidCB,'lon',ini);
lat_dim = netcdf.defDim(ncidCB,'lat',inj);
time_dim = netcdf.defDim(ncidCB,'time',nt);

vidlat = netcdf.defVar(ncidCB,'lat','double',lat_dim);
netcdf.putAtt(ncidCB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidCB,vidlat,'standard_name','lat');
netcdf.putAtt(ncidCB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidCB,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidCB,'lon','double',lon_dim);
netcdf.putAtt(ncidCB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidCB,vidlon,'standard_name','lon');
netcdf.putAtt(ncidCB,vidlon,'units','degrees_east');
netcdf.putAtt(ncidCB,vidlon,'axis','X');

vidtCB = netcdf.defVar(ncidCB,'time','double',time_dim);
netcdf.putAtt(ncidCB,vidtCB,'long_name','time');
netcdf.putAtt(ncidCB,vidtCB,'standard_name','time');
netcdf.putAtt(ncidCB,vidtCB,'calendar','365_day');
netcdf.putAtt(ncidCB,vidtCB,'axis','T');
netcdf.putAtt(ncidCB,vidtCB,'units','year');

vidbioCB = netcdf.defVar(ncidCB,'tcb','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidCB,vidbioCB,'long_name','total consumer biomass');
netcdf.putAtt(ncidCB,vidbioCB,'units','grams wet weight m-2' );
netcdf.defVarFill(ncidCB,vidbioCB,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidCB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidCB,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidCB,varid,'contact','C. Petrik');
netcdf.putAtt(ncidCB,varid,'institution','UC San Diego');
netcdf.putAtt(ncidCB,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidCB,varid,'includes benthos','yes');

netcdf.endDef(ncidCB);

netcdf.putVar(ncidCB,vidlat,ILAT);
netcdf.putVar(ncidCB,vidlon,ILON);
netcdf.putVar(ncidCB,vidbioCB,tcbI);
netcdf.putVar(ncidCB,vidtCB,yr);

netcdf.close(ncidCB);

%% tcp IPSL
ncidSB = netcdf.create(file_tcpI,'netcdf4');

lon_dim = netcdf.defDim(ncidSB,'lon',ini);
lat_dim = netcdf.defDim(ncidSB,'lat',inj);
time_dim = netcdf.defDim(ncidSB,'time',nt);

vidlat = netcdf.defVar(ncidSB,'lat','double',lat_dim);
netcdf.putAtt(ncidSB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'standard_name','lat');
netcdf.putAtt(ncidSB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidSB,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidSB,'lon','double',lon_dim);
netcdf.putAtt(ncidSB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidSB,vidlon,'standard_name','lon');
netcdf.putAtt(ncidSB,vidlon,'units','degrees_east');
netcdf.putAtt(ncidSB,vidlon,'axis','X');

vidtSB = netcdf.defVar(ncidSB,'time','double',time_dim);
netcdf.putAtt(ncidSB,vidtSB,'long_name','time');
netcdf.putAtt(ncidSB,vidtSB,'standard_name','time');
netcdf.putAtt(ncidSB,vidtSB,'units','year');
netcdf.putAtt(ncidSB,vidtSB,'calendar','365_day');
netcdf.putAtt(ncidSB,vidtSB,'axis','T');

vidbioSB = netcdf.defVar(ncidSB,'tpb','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSB,vidbioSB,'long_name','total consumer productivity');
netcdf.putAtt(ncidSB,vidbioSB,'units','grams wet weight m-2 y-1' );
netcdf.defVarFill(ncidSB,vidbioSB,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSB,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidSB,varid,'contact','C. Petrik');
netcdf.putAtt(ncidSB,varid,'institution','UC San Diego');
netcdf.putAtt(ncidSB,varid,'wet weight:C ratio','9:1');
netcdf.putAtt(ncidSB,varid,'includes benthos','no');

netcdf.endDef(ncidSB);

netcdf.putVar(ncidSB,vidlat,ILAT);
netcdf.putVar(ncidSB,vidlon,ILON);
netcdf.putVar(ncidSB,vidbioSB,tcpI);
netcdf.putVar(ncidSB,vidtSB,yr);

netcdf.close(ncidSB);



