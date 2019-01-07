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
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

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
cfileA = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE03_noCC_RE00100';
fpathA=['/Volumes/GFDL/NC/Matlab_new_size/' cfileA '/'];
load([fpathA 'Climatol_All_fish03.mat'],...
    'mp_my','lp_my','md_my','ld_my');
AallP=mp_my+lp_my;
AallD=md_my+ld_my;
AP=NaN*ones(ni,nj);
AD=NaN*ones(ni,nj);
AP(ID)=AallP;
AD(ID)=AallD;
clear mp_my lp_my md_my ld_my

cfileB = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathB=['/Volumes/GFDL/NC/Matlab_new_size/' cfileB '/'];
load([fpathB 'Climatol_All_fish03.mat'],...
    'mp_my','lp_my','md_my','ld_my');
BallP=mp_my+lp_my;
BallD=md_my+ld_my;
BP=NaN*ones(ni,nj);
BD=NaN*ones(ni,nj);
BP(ID)=BallP;
BD(ID)=BallD;
clear mp_my lp_my md_my ld_my

cfileC = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpathC=['/Volumes/GFDL/NC/Matlab_new_size/' cfileC '/'];
ppath = [pp cfileC '/'];
load([fpathC 'Means_bio_prod_fish_Climatol_All_fish03_',cfileC,'.mat'],...
    'mp_my','lp_my','md_my','ld_my');
CallP=mp_my+lp_my;
CallD=md_my+ld_my;
CP=NaN*ones(ni,nj);
CD=NaN*ones(ni,nj);
CP(ID)=CallP;
CD(ID)=CallD;
clear mp_my lp_my md_my ld_my

cfileD = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE10_noCC_RE00100';
fpathD=['/Volumes/GFDL/NC/Matlab_new_size/' cfileD '/'];
load([fpathD 'Climatol_All_fish03.mat'],...
    'mp_my','lp_my','md_my','ld_my');
DallP=mp_my+lp_my;
DallD=md_my+ld_my;
DP=NaN*ones(ni,nj);
DD=NaN*ones(ni,nj);
DP(ID)=DallP;
DD(ID)=DallD;
clear mp_my lp_my md_my ld_my

cfileE = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE13_noCC_RE00100';
fpathE=['/Volumes/GFDL/NC/Matlab_new_size/' cfileE '/'];
load([fpathE 'Climatol_All_fish03.mat'],...
    'mp_my','lp_my','md_my','ld_my');
EallP=mp_my+lp_my;
EallD=md_my+ld_my;
EP=NaN*ones(ni,nj);
ED=NaN*ones(ni,nj);
EP(ID)=EallP;
ED(ID)=EallD;
clear mp_my lp_my md_my ld_my

%% P/(P+D)
AfracPD = AP ./ (AP+AD);
BfracPD = BP ./ (BP+BD);
CfracPD = CP ./ (CP+CD);
DfracPD = DP ./ (DP+DD);
EfracPD = EP ./ (EP+ED);

APcatch = AP .* AREA_OCN * 365;
ADcatch = AD .* AREA_OCN * 365;
BPcatch = BP .* AREA_OCN * 365;
BDcatch = BD .* AREA_OCN * 365;
CPcatch = CP .* AREA_OCN * 365;
CDcatch = CD .* AREA_OCN * 365;
DPcatch = DP .* AREA_OCN * 365;
DDcatch = DD .* AREA_OCN * 365;
EPcatch = EP .* AREA_OCN * 365;
EDcatch = ED .* AREA_OCN * 365;

tlme = lme_mask_onedeg;
lme_A = NaN*ones(66,2);
lme_B = NaN*ones(66,2);
lme_C = NaN*ones(66,2);
lme_D = NaN*ones(66,2);
lme_E = NaN*ones(66,2);
lme_area = NaN*ones(66,1);
lme_lat = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    %total catch g
    lme_A(L,1) = nansum(APcatch(lid));
    lme_A(L,2) = nansum(ADcatch(lid));
    lme_B(L,1) = nansum(BPcatch(lid));
    lme_B(L,2) = nansum(BDcatch(lid));
    lme_C(L,1) = nansum(CPcatch(lid));
    lme_C(L,2) = nansum(CDcatch(lid));
    lme_D(L,1) = nanmean(DPcatch(lid));
    lme_D(L,2) = nanmean(DDcatch(lid));
    lme_E(L,1) = nanmean(EPcatch(lid));
    lme_E(L,2) = nanmean(EDcatch(lid));
    %total area of LME
    lme_area(L) = nansum(AREA_OCN(lid));
    %mean latt of LME
    lme_lat(L) = nanmean(geolat_t(lid));
end

lme_area_km2 = repmat(lme_area,1,2) * 1e-6;
% MT/km2
lme_A2 = lme_A * 1e-6 ./ lme_area_km2;
lme_B2 = lme_B * 1e-6 ./ lme_area_km2;
lme_C2 = lme_C * 1e-6 ./ lme_area_km2;
lme_D2 = lme_D * 1e-6 ./ lme_area_km2;
lme_E2 = lme_E * 1e-6 ./ lme_area_km2;

AlmePD = lme_A2(:,1) ./ sum(lme_A2,2);
BlmePD = lme_B2(:,1) ./ sum(lme_B2,2);
ClmePD = lme_C2(:,1) ./ sum(lme_C2,2);
DlmePD = lme_D2(:,1) ./ sum(lme_D2,2);
ElmePD = lme_E2(:,1) ./ sum(lme_E2,2);

grid_APD = NaN*ones(180,360);
grid_BPD = grid_APD;
grid_CPD = grid_APD;
grid_DPD = grid_APD;
grid_EPD = grid_APD;

for L=1:66
    lid = find(tlme==L);
    grid_APD(lid) = AlmePD(L,1);
    grid_BPD(lid) = BlmePD(L,1);
    grid_CPD(lid) = ClmePD(L,1);
    grid_DPD(lid) = DlmePD(L,1);
    grid_EPD(lid) = ElmePD(L,1);
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
print('-dpng',[ppath 'Climatol_' harv '_global_PvsD_comp_BE.png'])

%% 
figure(2)
%A
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,grid_APD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
text(-0.75,1.75,'\beta=0.025')
text(-2.75,1.75,'A')
%
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,grid_BPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
text(-0.75,1.75,'\beta=0.05')
text(-2.75,1.75,'B')
%
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,grid_CPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
text(-0.75,1.75,'\beta=0.075')
text(-2.75,1.75,'C')
%
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,grid_DPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
text(-0.75,1.75,'\beta=0.10')
text(-2.75,1.75,'D')
%
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,grid_EPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.4375 0.32 0.34 0.025],'orientation','horizontal','AxisLocation','in')
set(gcf,'renderer','painters')
text(-0.75,1.75,'\beta=0.125')
text(-2.75,1.75,'E')
stamp(cfileC)

% subplot('Position',[0.41 0.06 0.4 0.3])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,grid_FPD)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0 1]);
% set(gcf,'renderer','painters')
% text(-2.75,1.75,'\beta=0.15')
%stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_lme_PvsD_comp_BE_ms.png'])

%%
figure(3)
plot(AlmePD(:),lme_lat(:),'.k','MarkerSize',10); hold on;
plot(BlmePD(:),lme_lat(:),'+r','MarkerSize',10); hold on;
plot(ClmePD(:),lme_lat(:),'^b','MarkerSize',10); hold on;
plot(DlmePD(:),lme_lat(:),'og','MarkerSize',10); hold on;
plot(ElmePD(:),lme_lat(:),'*c','MarkerSize',10); hold on;
ylim([-90 90])

[slme_lat,ix] = sort(lme_lat);
rmat(:,1) = lme_lat(ix);
rmat(:,2) = AlmePD(ix);
rmat(:,3) = BlmePD(ix);
rmat(:,4) = ClmePD(ix);
rmat(:,5) = DlmePD(ix);
rmat(:,6) = ElmePD(ix);
rmat(:,7) = round(slme_lat,-1);

rtab = array2table(rmat,'VariableNames',{'LMElat','B025','B050','B075',...
    'B100','B125','RoundLat'});
writetable(rtab,'/Volumes/GFDL/NC/Matlab_new_size/bio_rates/LME_Pfrac_BEs.csv','Delimiter',',')


