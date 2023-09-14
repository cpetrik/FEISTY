% Time series of NPP, Z, and Det
% Area-integrated totals
% Test depth-weighted ratio

clear 
close all

Pdir = '/Volumes/petrik-lab/Feisty/GCM_Data/ESM2M_hist/';
cpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/grid_cobalt/';
gpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/cobalt_data/';
pp = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

%% GCM grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t','AREA_OCN'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);
AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

% Temp, Zoop, det, and npp 
load([gpath 'cobalt_temp_means.mat'],'ptemp_1yr_hist','btemp_1yr_hist',...
    'ptemp_1yr_fore','btemp_1yr_fore');
load([gpath 'cobalt_det_biom_means.mat'],'det_1yr_hist','det_1yr_fore');
load([gpath 'cobalt_npp_means.mat'],'npp_1yr_hist','npp_1yr_fore');
load([gpath 'cobalt_zoop_biom_means.mat'],...
    'mzprod_1yr_hist','lzprod_1yr_hist','mzprod_1yr_fore','lzprod_1yr_fore');

%% Zoop and det and npp 

%ESM2M in mmol N m-2 or mmol N m-2 d-1
% molN/m2/s --> gC/m2/d
% 106/16 mol C in 1 mol N
% 12.01 g C in 1 mol C
mzprod_1yr_hist = mzprod_1yr_hist * (106.0/16.0) * 12.01 * 60 * 60 * 24;
lzprod_1yr_hist = lzprod_1yr_hist * (106.0/16.0) * 12.01 * 60 * 60 * 24;
det_1yr_hist = det_1yr_hist * (106.0/16.0) * 12.01 * 60 * 60 * 24;
npp_1yr_hist = npp_1yr_hist * (106.0/16.0) * 12.01 * 60 * 60 * 24;

hmmz_prod = mzprod_1yr_hist;
hmlz_prod = lzprod_1yr_hist;
hmdet = det_1yr_hist;
hmnpp = npp_1yr_hist;
hmz = hmmz_prod + hmlz_prod;

mzprod_1yr_fore = mzprod_1yr_fore * (106.0/16.0) * 12.01 * 60 * 60 * 24;
lzprod_1yr_fore = lzprod_1yr_fore * (106.0/16.0) * 12.01 * 60 * 60 * 24;
det_1yr_fore = det_1yr_fore * (106.0/16.0) * 12.01 * 60 * 60 * 24;
npp_1yr_fore = npp_1yr_fore * (106.0/16.0) * 12.01 * 60 * 60 * 24;

fmmz_prod = mzprod_1yr_fore;
fmlz_prod = lzprod_1yr_fore;
fmdet = det_1yr_fore;
fmnpp = npp_1yr_fore;
fmz = fmmz_prod + fmlz_prod;

ptemp_1yr_hist = ptemp_1yr_hist - 273;
ptemp_1yr_fore = ptemp_1yr_fore - 273;

%% time
y1 = 1860+(1/12):1:2005;
y2 = 2005+(1/12):1:2100;

%% find N Atl & Artic region
ilat = find(grid(:,3)>=40);
ilon = find(grid(:,2)>=-100);
inatl = intersect(ilat,ilon);

IDnatl = ID(inatl);
grid_natl = grid(inatl,:);

%% put everything on grid then cut away

[ni,nj]=size(geolon_t);

Zz = NaN*ones(ni,nj);
Zz(IDnatl) = hmz(inatl,1);

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
surfm(geolat_t,geolon_t,Zz)
colormap('parula')
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% put everything on grid then cut away
rnan = find(all(isnan(Zz),2));
cnan = find(all(isnan(Zz),1));

rids = setdiff(1:ni,rnan);
cids = setdiff(1:nj,cnan);

tlat = geolat_t(rids,cids);
tlon = geolon_t(rids,cids);

%% put everything on grid then cut away
nyh=length(y1);

HTP = NaN*ones(ni,nj,nyh);
HTB = NaN*ones(ni,nj,nyh);
HDet = NaN*ones(ni,nj,nyh);
HNPP = NaN*ones(ni,nj,nyh);
HZM = NaN*ones(ni,nj,nyh);

for t=1:nyh
    Hp=NaN*ones(ni,nj);
    Hz=NaN*ones(ni,nj);
    Htp=NaN*ones(ni,nj);
    Htb=NaN*ones(ni,nj);
    Hdet=NaN*ones(ni,nj);

    Hp(IDnatl) = hmnpp(inatl,t);
    Hz(IDnatl) = hmz(inatl,t);
    Htp(IDnatl) = ptemp_1yr_hist(inatl,t);
    Htb(IDnatl) = btemp_1yr_hist(inatl,t);
    Hdet(IDnatl) = hmdet(inatl,t);

    HTP(:,:,t) = Htp;
    HTB(:,:,t) = Htb;
    HDet(:,:,t) = Hdet;
    HNPP(:,:,t) = Hp;
    HZM(:,:,t) = Hz;

end

%%
nyf=length(y2);

FTP = NaN*ones(ni,nj,nyf);
FTB = NaN*ones(ni,nj,nyf);
FDet = NaN*ones(ni,nj,nyf);
FNPP = NaN*ones(ni,nj,nyf);
FZM = NaN*ones(ni,nj,nyf);

for t=1:nyf
    Fp=NaN*ones(ni,nj);
    Fz=NaN*ones(ni,nj);
    Ftp=NaN*ones(ni,nj);
    Ftb=NaN*ones(ni,nj);
    Fdet=NaN*ones(ni,nj);

    Fp(IDnatl) = fmnpp(inatl,t);
    Fz(IDnatl) = fmz(inatl,t);
    Ftp(IDnatl) = ptemp_1yr_fore(inatl,t);
    Ftb(IDnatl) = btemp_1yr_fore(inatl,t);
    Fdet(IDnatl) = fmdet(inatl,t);

    FTP(:,:,t) = Htp;
    FTB(:,:,t) = Htb;
    FDet(:,:,t) = Hdet;
    FNPP(:,:,t) = Hp;
    FZM(:,:,t) = Hz;

end

%%
HTP = HTP(rids,cids,:);
HTB = HTB(rids,cids,:);
HDet = HDet(rids,cids,:);
HNPP = HNPP(rids,cids,:);
HZM = HZM(rids,cids,:);

FTP = FTP(rids,cids,:);
FTB = FTB(rids,cids,:);
FDet = FDet(rids,cids,:);
FNPP = FNPP(rids,cids,:);
FZM = FZM(rids,cids,:);

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
surfm(tlat,tlon,log10(squeeze(HDet(:,:,90))+eps))
colormap('parula')
colorbar
caxis([-4 0])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% save
save([gpath 'NAtl_ann_means_ESM2M_Hist_Forecast.mat'],'y1','y2',...
    'HTP','HTB','HDet','HNPP','HZM',...
    'FTP','FTB','FDet','FNPP','FZM',...
    'geolon_t','geolat_t','tlat','tlon',...
    'grid','grid_natl','inatl','ID','IDnatl');

spath = '/Volumes/petrik-lab/Feisty/GCM_Data/';
save([spath 'NAtl_ann_means_ESM2M_Hist_Forecast.mat'],'y1','y2',...
    'HTP','HTB','HDet','HNPP','HZM',...
    'FTP','FTB','FDet','FNPP','FZM',...
    'geolon_t','geolat_t','tlat','tlon',...
    'grid','grid_natl','inatl','ID','IDnatl');



