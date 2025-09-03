% Netcdf of LME mask for MOM6 COBALT2 sims
% prepared by Matthias at ISIMIP

clear 
close all

%%
Pdir = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';

ncid = netcdf.open([Pdir 'LMEs66_masks_1deg.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([Pdir 'cellarea_onedeg.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);


%% 
[LAT,LON] = meshgrid(lat,lon);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

%%
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,mask_NBS)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar;
clim([0 1])

%% 
% fuck this, I have to loop through all 66 masks to put on one grid???
Ldir = '/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/';

load([Ldir 'LMEs_ISIMIP_abbrev.mat'])

%%
[ni,nj] = size(cell_area);
lme_mask = nan*ones(ni,nj);

for i=1:length(LME)
    lname = char(varName(i));
    eval(['lid = find(' lname '(:)==1);']);
    lme_mask(lid) = LME(i);
end

%%
figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,lme_mask)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar;
clim([0 66])

%%
save([Ldir 'isimip_onedeg_lmemask_cellarea.mat'],'lme_mask','cell_area',...
    'LAT','LON','lat','lon');


