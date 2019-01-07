% Visualize output of POEM
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/Bio_rates/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
AREA_OCN = max(area,1);

%
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

close all

% plot info
[ni,nj]=size(lon);
geolat_t=lat;
geolon_t=lon;
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));


% colors
load('MyColormaps.mat')
cm9=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...  %b
    0 0 0];...      %black
    
cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm9);

%% Plots in space
cfileA = 'Dc_enc70-b200_cm20_m-b175-k06_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathA=['/Volumes/GFDL/CSV/Matlab_new_size/' cfileA '/'];
load([fpathA 'Clim_means_All_fish03.mat'],...
    'mf_my','mp_my','lp_my','md_my','ld_my');
AallF=mf_my;
AallP=mp_my+lp_my;
AallD=md_my+ld_my;
AF=NaN*ones(ni,nj);
AP=NaN*ones(ni,nj);
AD=NaN*ones(ni,nj);
AF(ID)=AallF;
AP(ID)=AallP;
AD(ID)=AallD;
clear mf_my mp_my lp_my md_my ld_my

cfileB = 'Dc_enc70-b200_cm20_m-b175-k08_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathB=['/Volumes/GFDL/CSV/Matlab_new_size/' cfileB '/'];
load([fpathB 'Clim_means_All_fish03.mat'],...
    'mf_my','mp_my','lp_my','md_my','ld_my');
BallF=mf_my;
BallP=mp_my+lp_my;
BallD=md_my+ld_my;
BF=NaN*ones(ni,nj);
BP=NaN*ones(ni,nj);
BD=NaN*ones(ni,nj);
BF(ID)=BallF;
BP(ID)=BallP;
BD(ID)=BallD;
clear mf_my mp_my lp_my md_my ld_my

cfileC = 'Dc_enc70-b200_cm20_m-b175-k10_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathC=['/Volumes/GFDL/CSV/Matlab_new_size/' cfileC '/'];
load([fpathC 'Clim_means_All_fish03.mat'],...
    'mf_my','mp_my','lp_my','md_my','ld_my');
CallF=mf_my;
CallP=mp_my+lp_my;
CallD=md_my+ld_my;
CF=NaN*ones(ni,nj);
CP=NaN*ones(ni,nj);
CD=NaN*ones(ni,nj);
CF(ID)=CallF;
CP(ID)=CallP;
CD(ID)=CallD;
clear mf_my mp_my lp_my md_my ld_my

cfileD = 'Dc_enc70-b200_cm20_m-b175-k12_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathD=['/Volumes/GFDL/CSV/Matlab_new_size/' cfileD '/'];
load([fpathD 'Clim_means_All_fish03.mat'],...
    'mf_my','mp_my','lp_my','md_my','ld_my');
DallF=mf_my;
DallP=mp_my+lp_my;
DallD=md_my+ld_my;
DF=NaN*ones(ni,nj);
DP=NaN*ones(ni,nj);
DD=NaN*ones(ni,nj);
DP(ID)=DallF;
DP(ID)=DallP;
DD(ID)=DallD;
clear mf_my mp_my lp_my md_my ld_my

%% P/(P+D)
AfracPD = AP ./ (AP+AD);
BfracPD = BP ./ (BP+BD);
CfracPD = CP ./ (CP+CD);
DfracPD = DP ./ (DP+DD);

AFcatch = AF .* AREA_OCN * 365;
APcatch = AP .* AREA_OCN * 365;
ADcatch = AD .* AREA_OCN * 365;
BFcatch = BF .* AREA_OCN * 365;
BPcatch = BP .* AREA_OCN * 365;
BDcatch = BD .* AREA_OCN * 365;
CFcatch = CF .* AREA_OCN * 365;
CPcatch = CP .* AREA_OCN * 365;
CDcatch = CD .* AREA_OCN * 365;
DFcatch = DF .* AREA_OCN * 365;
DPcatch = DP .* AREA_OCN * 365;
DDcatch = DD .* AREA_OCN * 365;

tlme = lme_mask_onedeg;
lme_A = NaN*ones(66,3);
lme_B = NaN*ones(66,3);
lme_C = NaN*ones(66,3);
lme_D = NaN*ones(66,3);
lme_area = NaN*ones(66,1);
lme_lat = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    %total catch g
    lme_A(L,1) = nansum(AFcatch(lid));
    lme_A(L,2) = nansum(APcatch(lid));
    lme_A(L,3) = nansum(ADcatch(lid));
    lme_B(L,1) = nansum(BFcatch(lid));
    lme_B(L,2) = nansum(BPcatch(lid));
    lme_B(L,3) = nansum(BDcatch(lid));
    lme_C(L,1) = nansum(CFcatch(lid));
    lme_C(L,2) = nansum(CPcatch(lid));
    lme_C(L,3) = nansum(CDcatch(lid));
    lme_D(L,1) = nanmean(DFcatch(lid));
    lme_D(L,2) = nanmean(DPcatch(lid));
    lme_D(L,3) = nanmean(DDcatch(lid));
    %total area of LME
    lme_area(L) = nansum(AREA_OCN(lid));
    %mean latt of LME
    lme_lat(L) = nanmean(geolat_t(lid));
end

lme_area_km2 = repmat(lme_area,1,3) * 1e-6;
% MT/km2
lme_A2 = lme_A * 1e-6 ./ lme_area_km2;
lme_B2 = lme_B * 1e-6 ./ lme_area_km2;
lme_C2 = lme_C * 1e-6 ./ lme_area_km2;
lme_D2 = lme_D * 1e-6 ./ lme_area_km2;

AlmePD = lme_A2(:,2) ./ (lme_A2(:,2) + lme_A2(:,3));
BlmePD = lme_B2(:,2) ./ (lme_B2(:,2) + lme_B2(:,3));
ClmePD = lme_C2(:,2) ./ (lme_C2(:,2) + lme_C2(:,3));
DlmePD = lme_D2(:,2) ./ (lme_D2(:,2) + lme_D2(:,3));

grid_APD = NaN*ones(180,360);
grid_BPD = grid_APD;
grid_CPD = grid_APD;
grid_DPD = grid_APD;

for L=1:66
    lid = find(tlme==L);
    grid_APD(lid) = AlmePD(L,1);
    grid_BPD(lid) = BlmePD(L,1);
    grid_CPD(lid) = ClmePD(L,1);
    grid_DPD(lid) = DlmePD(L,1);
end

%% 
figure(1)
%A
subplot('Position',[0 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,AfracPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'A')

%B
subplot('Position',[0.5 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,BfracPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'B')

%C
subplot('Position',[0 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,CfracPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'C')

%D
subplot('Position',[0.5 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,DfracPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.25 0.075 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
text(-2.75,1.75,'D')
%stamp([harv '_' cfile])
print('-dpng',[pp 'Climatol_' harv '_global_PvsD_comp_kt_6-12.png'])

%% 
figure(2)
%A
subplot('Position',[0 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,grid_APD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'A')

%B
subplot('Position',[0.5 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,grid_BPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'B')

%C
subplot('Position',[0 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,grid_CPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'C')

%D
subplot('Position',[0.5 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,grid_DPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.25 0.075 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
text(-2.75,1.75,'D')
%stamp([harv '_' cfile])
print('-dpng',[pp 'Climatol_' harv '_lme_PvsD_comp_kt_6-12.png'])

%% 

