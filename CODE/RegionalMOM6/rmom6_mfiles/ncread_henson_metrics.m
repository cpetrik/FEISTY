% S Henson 2012 results

clear 
close all

%%
epath = '/Volumes/petrik-lab/Feisty/Obs_data/Export_stat_products/S_Henson/Henson2012_Figure4_ncs/';

ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/OCB_BGC/export_flux/Henson2012/';

ncdisp([epath 'export_Henson2012.nc'])

%%
% EPCHn
% Size:       360x200
% Dimensions: lon,lat
% Datatype:   double
% Attributes:
EP_units      = 'gC m-2 yr-1';
% _FillValue = -999
% lat
% Size:       200x1
% Dimensions: lat
% Datatype:   double
% lon
% Size:       360x1
% Dimensions: lon
% Datatype:   double

%% export flux  (at 100 m?)
ncid = netcdf.open([epath 'export_Henson2012.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

EPCHn(EPCHn<0)=nan;

%% map info
[geolat,geolon] = meshgrid(lat,lon);

[ni,nj]=size(geolon);
geolon = double(geolon);
geolat = double(geolat);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;  

%% maps

test = squeeze(EPCHn(:,:,1));

figure; pcolor(geolat); shading flat; colorbar
figure; pcolor(geolon); shading flat; colorbar
figure; pcolor(test); shading flat; colorbar

%%
close all
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,EPCHn)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%clim([-3 0]);
hcb = colorbar('h');
title('mean export flux (g/m2/yr)')
stamp('')
print('-dpng',[ppath 'export_Henson2012_gCyr.png'])

%%
figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(EPCHn))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%clim([-3 0]);
hcb = colorbar('h');
title('log10 mean export flux (g/m2/yr)')
stamp('')
print('-dpng',[ppath 'export_log10_Henson2012_gCyr.png'])

%% same units as Pippa Edwards
EmgCd = EPCHn*1e3/365;

figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(EmgCd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 3]);
hcb = colorbar('h');
title('log10 mean export flux (mg/m2/day)')
stamp('')
print('-dpng',[ppath 'export_log10_Henson2012_mgCd.png'])

%% ------------------------Martin b-------------------------------
ncdisp([epath 'Martinb_Henson2012.nc'])

ncid = netcdf.open([epath 'Martinb_Henson2012.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Martin_b(Martin_b<-900)=nan;

%% ------------------------PE effic-------------------------------
ncdisp([epath 'PEeff_Henson2012.nc'])

ncid = netcdf.open([epath 'PEeff_Henson2012.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

pe(pe<0)=nan;

%% ------------------------poc2000-------------------------------
ncdisp([epath 'poc2000_Henson2012.nc'])
% pocC_L
% Size:       360x200
% Dimensions: lon,lat
% Datatype:   double
% Attributes:
pocC_L_units      = 'gC m-2 yr-1';
% _FillValue = -999

ncid = netcdf.open([epath 'poc2000_Henson2012.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

pocC_L(pocC_L<0)=nan;

%% ------------------------Teff-------------------------------
ncdisp([epath 'Teff_Henson2012.nc'])

ncid = netcdf.open([epath 'Teff_Henson2012.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

teffCHL(teffCHL<0)=nan;

%% save all 
save([epath 'Henson2012_all_vars.mat'],'EPCHn','EP_units',...
    'Martin_b','pe','pocC_L','pocC_L_units','teffCHL',...
    'lat','lon')