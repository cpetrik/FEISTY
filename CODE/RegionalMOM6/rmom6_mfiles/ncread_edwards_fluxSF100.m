% Pippa Edwards POC flux at seafloor

clear 
close all

%%
epath = '/Volumes/petrik-lab/Feisty/Globus_RO/export_obs_products/NWA_regrid/';

ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/OCB_BGC/export_flux/Edwards/';

ncdisp([epath 'nwa_estpocbayes_poc_SF100_mapped.nc'])

%%
meanlogpoc_units       = 'ln(mg/m2/day)';
meanlogpoc_descrip = 'Log-transformed mean poc generated from mu and sigma';
meanpoc_units       = 'mg/m2/day';
meanpoc_descrip = 'Mean poc generated from mu and sigma then taken out of log-space';
depth_units       = 'm';
depth_descrip = 'Depth at which flux is calculated at';
date_units    = 'days since 1997-09-01';
date_calendar = 'proleptic_gregorian';

%% POC flux seafloor
ncid = netcdf.open([epath 'nwa_estpocbayes_poc_SF100_mapped.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%% map info
[ni,nj]=size(lon);
geolon = double(lon);
geolat = double(lat);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;  

%% maps
yr=(double(date)/365)+1997+(9/12);
% mean over whole time
mpoc = mean(meanlogpoc,3,'omitnan');

figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,mpoc)
cmocean('matter')
%h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([2 5]);
hcb = colorbar('h');
title('log10 mean POC flux at seafloor (mg/m2/day) 1997-2022')
stamp('')
print('-dpng',[ppath 'nwa_meanlogpoc_SF100_mean.png'])

%% save
save([epath 'nwa_estpocbayes_poc_SF100_mapped.mat'],'date','date_units','date_calendar','yr',...
    'depth','depth_descrip','meanlogpoc','meanlogpoc_descrip','meanlogpoc_units',...
    'meanpoc','meanpoc_descrip','meanpoc_units','lat','lon')