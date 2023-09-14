% FEISTY output at all locations
% Annual means
% Just N Atl region

clear 
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/Historic_ESM2M/'];
spath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/NAtl/'];

%% SP
ncid = netcdf.open([fpath 'Historic_' harv '_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(biomass);
%%
SP.bio = biomass;
clear biomass prod

% SF
ncid = netcdf.open([fpath 'Historic_' harv '_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
clear biomass prod

% SD
ncid = netcdf.open([fpath 'Historic_' harv '_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
clear biomass prod

% MP
ncid = netcdf.open([fpath 'Historic_' harv '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
clear biomass yield

% MF
ncid = netcdf.open([fpath 'Historic_' harv '_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
clear biomass yield

% MD
ncid = netcdf.open([fpath 'Historic_' harv '_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
clear biomass yield

% LP
ncid = netcdf.open([fpath 'Historic_' harv '_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
clear biomass yield

% LD
ncid = netcdf.open([fpath 'Historic_' harv '_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
clear biomass yield

% Benthic material
ncid = netcdf.open([fpath 'Historic_' harv '_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
clear biomass

%% Group fn types

Fbio = SF.bio + MF.bio;
Pbio = SP.bio + MP.bio + LP.bio;
Dbio = SD.bio + MD.bio + LD.bio;

clear SF SP SD MF MP MD LP LD

%% Take means 

% Every 1 year
st=1:12:length(time);
en=12:12:length(time);

for n=1:length(st)
    f_mean(:,n)=mean(Fbio(:,st(n):en(n)),2,'omitnan');
    p_mean(:,n)=mean(Pbio(:,st(n):en(n)),2,'omitnan');
    d_mean(:,n)=mean(Dbio(:,st(n):en(n)),2,'omitnan');
    b_mean(:,n)=mean(Bent.bio(:,st(n):en(n)),2,'omitnan');
    
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
load coastlines;                     %decent looking coastlines

% iids = [206:295];
% jids = [150:177];

% plot info
plotminlat=40; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-100;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%
figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,geolon_t)
colormap('parula')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%%
ilat = find(grid(:,3)>=40);
ilon = find(grid(:,2)>=-100);
inatl = intersect(ilat,ilon);

IDnatl = ID(inatl);
grid_natl = grid(inatl,:);

%% check
[ni,nj]=size(geolon_t);

Zf=NaN*ones(ni,nj);

Zf(IDnatl) = f_mean(inatl,100);
%Zf(ID) = f_mean(:,1);

%%
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Zf)
colormap('parula')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%%
f_mean=f_mean(inatl,:);
p_mean=p_mean(inatl,:);
d_mean=d_mean(inatl,:);
b_mean=b_mean(inatl,:);

%% put everything on grid then cut away
rnan = find(all(isnan(Zf),2));
cnan = find(all(isnan(Zf),1));

rids = setdiff(1:ni,rnan);
cids = setdiff(1:nj,cnan);

test = Zf(rids,cids);
tlat = geolat_t(rids,cids);
tlon = geolon_t(rids,cids);

figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(tlat,tlon,test)
colormap('parula')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%%
ny=time(end)/12;

Fall = NaN*ones(ni,nj,ny);
Pall = NaN*ones(ni,nj,ny);
Dall = NaN*ones(ni,nj,ny);
Ball = NaN*ones(ni,nj,ny);

for t=1:ny
    Zf=NaN*ones(ni,nj);
    Zp=NaN*ones(ni,nj);
    Zd=NaN*ones(ni,nj);
    Zb=NaN*ones(ni,nj);

    Zf(IDnatl) = f_mean(:,t);
    Zp(IDnatl) = p_mean(:,t);
    Zd(IDnatl) = d_mean(:,t);
    Zb(IDnatl) = b_mean(:,t);

    Fall(:,:,t) = Zf;
    Pall(:,:,t) = Zp;
    Dall(:,:,t) = Zd;
    Ball(:,:,t) = Zb;

end

%%
Fall = Fall(rids,cids,:);
Pall = Pall(rids,cids,:);
Dall = Dall(rids,cids,:);
Ball = Ball(rids,cids,:);


%%
save([spath 'NAtl_Means_Historic_' harv '_' cfile '.mat'],'time',...
    'f_mean','p_mean','d_mean','b_mean','geolon_t','geolat_t',...
    'Fall','Pall','Dall','Ball','tlat','tlon',...
    'grid','grid_natl','inatl','ID','IDnatl');








