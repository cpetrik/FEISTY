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
load([cpath 'LME_clim_temp.mat']);
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

% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

x=-8:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% SAUP data
spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);
l10sD=log10(Dlme_mcatch10+eps);
sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

%% Plots in space
cfileA = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D025_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathA=['/Volumes/GFDL/NC/Matlab_new_size/' cfileA '/'];
load([fpathA 'Means_bio_prod_fish_Climatol_All_fish03_',cfileA,'.mat'],...
    'mp_my','lp_my','md_my','ld_my');
AallP=mp_my+lp_my;
AallD=md_my+ld_my;
AP=NaN*ones(ni,nj);
AD=NaN*ones(ni,nj);
AP(ID)=AallP;
AD(ID)=AallD;
clear mp_my lp_my md_my ld_my

cfileB = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D050_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathB=['/Volumes/GFDL/NC/Matlab_new_size/' cfileB '/'];
load([fpathB 'Means_bio_prod_fish_Climatol_All_fish03_',cfileB,'.mat'],...
    'mp_my','lp_my','md_my','ld_my');
BallP=mp_my+lp_my;
BallD=md_my+ld_my;
BP=NaN*ones(ni,nj);
BD=NaN*ones(ni,nj);
BP(ID)=BallP;
BD(ID)=BallD;
clear mp_my lp_my md_my ld_my

cfileC = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathC=['/Volumes/GFDL/NC/Matlab_new_size/' cfileC '/'];
load([fpathC 'Means_bio_prod_fish_Climatol_All_fish03_',cfileC,'.mat'],...
    'mp_my','lp_my','md_my','ld_my');
CallP=mp_my+lp_my;
CallD=md_my+ld_my;
CP=NaN*ones(ni,nj);
CD=NaN*ones(ni,nj);
CP(ID)=CallP;
CD(ID)=CallD;
clear mp_my lp_my md_my ld_my

cfileD = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D100_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathD=['/Volumes/GFDL/NC/Matlab_new_size/' cfileD '/'];
load([fpathD 'Means_bio_prod_fish_Climatol_All_fish03_',cfileD,'.mat'],...
    'mp_my','lp_my','md_my','ld_my');
DallP=mp_my+lp_my;
DallD=md_my+ld_my;
DP=NaN*ones(ni,nj);
DD=NaN*ones(ni,nj);
DP(ID)=DallP;
DD(ID)=DallD;
clear mp_my lp_my md_my ld_my

%% P/(P+D)
AfracPD = AP ./ (AP+AD);
BfracPD = BP ./ (BP+BD);
CfracPD = CP ./ (CP+CD);
DfracPD = DP ./ (DP+DD);

APcatch = AP .* AREA_OCN * 365;
ADcatch = AD .* AREA_OCN * 365;
BPcatch = BP .* AREA_OCN * 365;
BDcatch = BD .* AREA_OCN * 365;
CPcatch = CP .* AREA_OCN * 365;
CDcatch = CD .* AREA_OCN * 365;
DPcatch = DP .* AREA_OCN * 365;
DDcatch = DD .* AREA_OCN * 365;

tlme = lme_mask_onedeg;
lme_A = NaN*ones(66,2);
lme_B = NaN*ones(66,2);
lme_C = NaN*ones(66,2);
lme_D = NaN*ones(66,2);
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
    lme_D(L,1) = nansum(DPcatch(lid));
    lme_D(L,2) = nansum(DDcatch(lid));
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

AlmePD = lme_A2(:,1) ./ sum(lme_A2,2);
BlmePD = lme_B2(:,1) ./ sum(lme_B2,2);
ClmePD = lme_C2(:,1) ./ sum(lme_C2,2);
DlmePD = lme_D2(:,1) ./ sum(lme_D2,2);

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

Al10pD=log10(lme_A2(:,2));
Bl10pD=log10(lme_B2(:,2));
Cl10pD=log10(lme_C2(:,2));
Dl10pD=log10(lme_D2(:,2));

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
text(-2.75,1.75,'D=0.25')

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
text(-2.75,1.75,'D=0.50')

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
text(-2.75,1.75,'D=0.75')

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
text(-2.75,1.75,'D=1.0')
%stamp([harv '_' cfile])
print('-dpng',[pp 'Climatol_' harv '_global_PvsD_comp_Dpref.png'])

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
text(-2.75,1.75,'D=0.25')

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
text(-2.75,1.75,'D=0.50')

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
text(-2.75,1.75,'D=0.75')

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
text(-2.75,1.75,'D=1.0')
stamp('')
print('-dpng',[pp 'Climatol_' harv '_lme_PvsD_comp_Dpref.png'])

%% Do in R as a box-whisker plot
figure(3)
plot(AlmePD(:),lme_lat(:),'.k','MarkerSize',10); hold on;
plot(BlmePD(:),lme_lat(:),'+r','MarkerSize',10); hold on;
plot(ClmePD(:),lme_lat(:),'^b','MarkerSize',10); hold on;
plot(DlmePD(:),lme_lat(:),'og','MarkerSize',10); hold on;
ylim([-90 90])
xlim([0 1])
ylabel('Latitude')
title('P/(P+D)')
legend('.25','.5','.75','1')
stamp('')
print('-dpng',[pp 'Climatol_' harv '_latitude_PvsD_comp_Dpref.png'])

[slme_lat,ix] = sort(lme_lat);
rmat(:,1) = lme_lat(ix);
rmat(:,2) = AlmePD(ix);
rmat(:,3) = BlmePD(ix);
rmat(:,4) = ClmePD(ix);
rmat(:,5) = DlmePD(ix);
rmat(:,6) = round(slme_lat,-1);

rtab = array2table(rmat,'VariableNames',{'LMElat','D25','D50','D75','D100'});
writetable(rtab,'/Volumes/GFDL/NC/Matlab_new_size/bio_rates/LME_Pfrac_Dprefs.csv','Delimiter',',')

%%
figure(33)
subplot(2,2,1)
boxplot(AlmePD,rmat(:,6))
%xlim([-90 90])
ylim([0 1])
xlabel('Latitude')
title('D=0.25')

subplot(2,2,2)
boxplot(BlmePD,rmat(:,6))
%xlim([-90 90])
ylim([0 1])
xlabel('Latitude')
title('D=0.5')

subplot(2,2,3)
boxplot(ClmePD,rmat(:,6))
%xlim([-90 90])
ylim([0 1])
xlabel('Latitude')
title('D=0.75')

subplot(2,2,4)
boxplot(DlmePD,rmat(:,6))
%xlim([-90 90])
ylim([0 1])
xlabel('Latitude')
title('D=1.0')
stamp('')
print('-dpng',[pp 'Climatol_' harv '_latitude_box_PvsD_comp_Dpref.png'])

%% SAU corrs
rA=corr(l10sD(keep),Al10pD(keep));
rB=corr(l10sD(keep),Bl10pD(keep));
rC=corr(l10sD(keep),Cl10pD(keep));
rD=corr(l10sD(keep),Dl10pD(keep));
rPA=corr(sFracPD(keep),AlmePD(keep));
rPB=corr(sFracPD(keep),BlmePD(keep));
rPC=corr(sFracPD(keep),ClmePD(keep));
rPD=corr(sFracPD(keep),DlmePD(keep));

figure(4);
subplot(2,2,1)
plot(x,x,'--k');hold on;
for j=1:length(keep)
    lme=keep(j);
    plot(l10sD(lme),Al10pD(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
end
text(-1.75,1.7,['r = ' sprintf('%2.2f',rA)])
axis([-2 2 -2 2])
xlabel('SAU')
ylabel('POEM')
title('D=0.25')

subplot(2,2,2)
plot(x,x,'--k');hold on;
for j=1:length(keep)
    lme=keep(j);
    plot(l10sD(lme),Bl10pD(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
end
text(-1.75,1.7,['r = ' sprintf('%2.2f',rB)])
axis([-2 2 -2 2])
xlabel('SAU')
ylabel('POEM')
title('D=0.50')

subplot(2,2,3)
plot(x,x,'--k');hold on;
for j=1:length(keep)
    lme=keep(j);
    plot(l10sD(lme),Cl10pD(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
end
text(-1.75,1.7,['r = ' sprintf('%2.2f',rC)])
axis([-2 2 -2 2])
xlabel('SAU')
ylabel('POEM')
title('D=0.75')

subplot(2,2,4)
plot(x,x,'--k');hold on;
for j=1:length(keep)
    lme=keep(j);
    plot(l10sD(lme),Dl10pD(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
end
text(-1.75,1.7,['r = ' sprintf('%2.2f',rD)])
axis([-2 2 -2 2])
xlabel('SAU')
ylabel('POEM')
title('D=1.0')
stamp('')
print('-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_Dpref_tests_scatterD.png'])


figure(5);
subplot(2,2,1)
plot(x,x,'--k');hold on;
for j=1:length(keep)
    lme=keep(j);
    plot(sFracPD(lme),AlmePD(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.8,0.35,['r = ' sprintf('%2.2f',rPA)])
axis([0 1 0 1])
xlabel('SAU')
ylabel('POEM')
title('D=0.25')

subplot(2,2,2)
plot(x,x,'--k');hold on;
for j=1:length(keep)
    lme=keep(j);
    plot(sFracPD(lme),BlmePD(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.8,0.35,['r = ' sprintf('%2.2f',rPB)])
axis([0 1 0 1])
xlabel('SAU')
ylabel('POEM')
title('D=0.50')

subplot(2,2,3)
plot(x,x,'--k');hold on;
for j=1:length(keep)
    lme=keep(j);
    plot(sFracPD(lme),ClmePD(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.8,0.35,['r = ' sprintf('%2.2f',rPC)])
axis([0 1 0 1])
xlabel('SAU')
ylabel('POEM')
title('D=0.75')

subplot(2,2,4)
plot(x,x,'--k');hold on;
for j=1:length(keep)
    lme=keep(j);
    plot(sFracPD(lme),DlmePD(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.8,0.35,['r = ' sprintf('%2.2f',rPD)])
axis([0 1 0 1])
xlabel('SAU')
ylabel('POEM')
title('D=1.0')
stamp('')
print('-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_Dpref_tests_scatterPD.png'])


