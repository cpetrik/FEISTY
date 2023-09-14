% FEISTY output at all locations

clear 
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/Forecast_RCP85_ESM2M/'];
spath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/NAtl/'];

%% SP
ncid = netcdf.open([fpath 'Forecast_',harv,'_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(prod);

%
SP.prod = prod;

clear prod 

%% SF
ncid = netcdf.open([fpath 'Forecast_',harv,'_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.prod = prod(:,1:nt);

clear prod 

% SD
ncid = netcdf.open([fpath 'Forecast_',harv,'_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.prod = prod;

clear prod 

% MP
ncid = netcdf.open([fpath 'Forecast_',harv,'_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.prod = prod;

clear prod 

% MF
ncid = netcdf.open([fpath 'Forecast_',harv,'_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.prod = prod;

clear prod 

% MD
ncid = netcdf.open([fpath 'Forecast_',harv,'_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.prod = prod;

clear prod 
% LP
ncid = netcdf.open([fpath 'Forecast_',harv,'_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.prod = prod;

clear prod 

% LD
ncid = netcdf.open([fpath 'Forecast_',harv,'_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.prod = prod;

clear prod 
% Benthic material
ncid = netcdf.open([fpath 'Forecast_',harv,'_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);


%% Group fn types

Fprod = SF.prod + MF.prod;
Pprod = SP.prod + MP.prod + LP.prod;
Dprod = SD.prod + MD.prod + LD.prod;

clear SF SP SD MF MP MD LP LD

%% Take means 

% Every 1 year
st=1:12:length(time);
en=12:12:length(time);

for n=1:length(st)
    f_mean(:,n)=mean(Fprod(:,st(n):en(n)),2,'omitnan');
    p_mean(:,n)=mean(Pprod(:,st(n):en(n)),2,'omitnan');
    d_mean(:,n)=mean(Dprod(:,st(n):en(n)),2,'omitnan');
    
end

%% find N Atl & Artic region
cpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/grid_cobalt/';
pp = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

ppath = [pp cfile '/NAtl/'];

load('/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

%%
ilat = find(grid(:,3)>=40);
ilon = find(grid(:,2)>=-100);
inatl = intersect(ilat,ilon);

IDnatl = ID(inatl);
grid_natl = grid(inatl,:);

%%
% f_mean=f_mean(inatl,:);
% p_mean=p_mean(inatl,:);
% d_mean=d_mean(inatl,:);
% b_mean=b_mean(inatl,:);

%% put everything on grid then cut away

[ni,nj]=size(geolon_t);

Zf = NaN*ones(ni,nj);
Zf(IDnatl) = f_mean(inatl,1);

%%
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;                   

figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Zf)
colormap('parula')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% put everything on grid then cut away
rnan = find(all(isnan(Zf),2));
cnan = find(all(isnan(Zf),1));

rids = setdiff(1:ni,rnan);
cids = setdiff(1:nj,cnan);

tlat = geolat_t(rids,cids);
tlon = geolon_t(rids,cids);

%% put everything on grid then cut away
ny=time(end)/12;

Fall = NaN*ones(ni,nj,ny);
Pall = NaN*ones(ni,nj,ny);
Dall = NaN*ones(ni,nj,ny);

for t=1:ny
    Zf=NaN*ones(ni,nj);
    Zp=NaN*ones(ni,nj);
    Zd=NaN*ones(ni,nj);

    Zf(IDnatl) = f_mean(inatl,t);
    Zp(IDnatl) = p_mean(inatl,t);
    Zd(IDnatl) = d_mean(inatl,t);

    Fall(:,:,t) = Zf;
    Pall(:,:,t) = Zp;
    Dall(:,:,t) = Zd;

end

%%
Fall = Fall(rids,cids,:);
Pall = Pall(rids,cids,:);
Dall = Dall(rids,cids,:);

%% check
% iids = [206:295];
% jids = [150:177];

% plot info
plotminlat=40; 
plotmaxlat=90;
plotminlon=-100;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(tlat,tlon,log10(squeeze(Dall(:,:,90))+eps))
colormap('parula')
colorbar
caxis([-4 0])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% save
save([spath 'NAtl_Means_prod_Forecast_' harv '_' cfile '.mat'],'time',...
    'f_mean','p_mean','d_mean','geolon_t','geolat_t',...
    'Fall','Pall','Dall','tlat','tlon',...
    'grid','grid_natl','inatl','ID','IDnatl');





