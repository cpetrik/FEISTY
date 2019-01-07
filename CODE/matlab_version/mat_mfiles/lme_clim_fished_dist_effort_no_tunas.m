% Estimate POEM harvest with different efforts scaled from land distance
% Climatology
% 150 years
% Saved as mat files
% Fish within certain distance from shore

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([cpath 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'esm26_1deg_shore_dist.mat']);
load([cpath 'LME_clim_temp_zoop_det.mat']);

AREA_OCN = max(area,1);
tlme = lme_mask_onedeg;
% "lfile" never changes, has lme areas
lfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
lpath = ['/Volumes/GFDL/NC/Matlab_new_size/' lfile '/'];
load([lpath 'LME_clim_fished03_' lfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

%SAUP
load([spath 'SAUP_LME_Catch_top10_Stock_no_tunas_bathy.mat']);
l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

%POEM
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
harv = 'All_fish03';
fpath = [dp cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);

%% plot info
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

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')
%Assign a color to each LME based on temp
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

x=-6:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% Plots on global grid
Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);

Cmf(ID)=mf_my;
Cmp(ID)=mp_my;
Cmd(ID)=md_my;
Clp(ID)=lp_my;
Cld(ID)=ld_my;

%% Effort scalings
eff = min_dist;
eff(eff>600) = 0;
eff(eff>0) = 1;
%Make LMEs 10 (Hawaii), 65 (Aleutian Is) = 1
lid10 = find(tlme==10);
lid65 = find(tlme==65);
eff(lid10) = 1;
eff(lid65) = 1;
cfile2 = '_effort_600km_no_tunas_bathy';
tit2 = 'Effort 600km';

% Distance
figure(10)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(min_dist))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2 2.7]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Distance to land')
stamp('')

% Effort 
figure(11)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,eff)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title(tit2)
stamp('')
% print('-dpng',[pp harv,cfile2,'_LME.png'])


%% apply effort scaling
%g/m2/d --> total g
Amf_mcatch = Cmf .* AREA_OCN * 365 .* eff; %mean fish catch per yr
Amp_mcatch = Cmp .* AREA_OCN * 365 .* eff;
Amd_mcatch = Cmd .* AREA_OCN * 365 .* eff;
Alp_mcatch = Clp .* AREA_OCN * 365 .* eff;
Ald_mcatch = Cld .* AREA_OCN * 365 .* eff;

%% Calc LMEs
lme_mcatch = NaN*ones(66,5);

for L=1:66
    lid = find(tlme==L);
    %total catch g
    lme_mcatch(L,1) = nansum(Amf_mcatch(lid));
    lme_mcatch(L,2) = nansum(Amp_mcatch(lid));
    lme_mcatch(L,3) = nansum(Amd_mcatch(lid));
    lme_mcatch(L,4) = nansum(Alp_mcatch(lid));
    lme_mcatch(L,5) = nansum(Ald_mcatch(lid));
end
% in MT
plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;

% MT/km2
plme_mcatch = plme_mcatch ./ lme_area_km2;
plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;
pFracPD = plme_Pmcatch ./ (plme_Pmcatch + plme_Dmcatch);

l10p=log10(plme_mcatch+eps);
l10pF=log10(plme_Fmcatch+eps);
l10pP=log10(plme_Pmcatch+eps);
l10pD=log10(plme_Dmcatch+eps);

% Catch in LMEs vs open ocean
%total g, not corrected by area
tot_catch = (nansum(Amf_mcatch(:)) + nansum(Amp_mcatch(:)) + nansum(Amd_mcatch(:))...
    + nansum(Alp_mcatch(:)) + nansum(Ald_mcatch(:))) * 1e-6;
frac_lme = nansum(nansum(lme_mcatch,2) * 1e-6) / tot_catch

%% Stats
%r
rall=corr(l10s(keep),l10p(keep));
rF=corr(l10sF(keep),l10pF(keep));
rP=corr(l10sP(keep),l10pP(keep));
rD=corr(l10sD(keep),l10pD(keep));
rPD=corr(sFracPD(keep),pFracPD(keep));

%root mean square error
o=l10s(keep);
p=l10p(keep);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=l10sF(keep);
p=l10pF(keep);
n = length(o);
num=nansum((p-o).^2);
rmseF = sqrt(num/n);

o=l10sP(keep);
p=l10pP(keep);
n = length(o);
num=nansum((p-o).^2);
rmseP = sqrt(num/n);

o=l10sD(keep);
p=l10pD(keep);
n = length(o);
num=nansum((p-o).^2);
rmseD = sqrt(num/n);

o=sFracPD(keep);
p=pFracPD(keep);
n = length(o);
num=nansum((p-o).^2);
rmsePD = sqrt(num/n);

%Fmed
Fall=10^(median(l10s(keep)-l10p(keep)));
FF=10^(median(l10sF(keep)-l10pF(keep)));
FP=10^(median(l10sP(keep)-l10pP(keep)));
FD=10^(median(l10sD(keep)-l10pD(keep)));
FPD=10^(median(sFracPD(keep)-pFracPD(keep)));

% Table
fish_stat(1,1) = rall;
fish_stat(2,1) = rF;
fish_stat(3,1) = rP;
fish_stat(4,1) = rD;
fish_stat(1,2) = rmse;
fish_stat(2,2) = rmseF;
fish_stat(3,2) = rmseP;
fish_stat(4,2) = rmseD;
fish_stat(1,3) = Fall;
fish_stat(2,3) = FF;
fish_stat(3,3) = FP;
fish_stat(4,3) = FD;
fish_stat(5,1) = rPD;
fish_stat(5,2) = rmsePD;
fish_stat(5,3) = FPD;

%% save
Fstat = array2table(fish_stat,'VariableNames',{'r','RMSE','Fmed'},...
    'RowNames',{'All','F','P','D','FracP'});
writetable(Fstat,[fpath 'LME_SAUP_stats' cfile2 '_' cfile '.csv'],'Delimiter',',',...
    'WriteRowNames',true)
save([fpath 'LME_SAUP_stats' cfile2 '_' cfile '.mat'],'fish_stat')

save([fpath 'LME_clim_',harv,cfile2,'.mat'],'lme_mcatch');

%% Figures
clme_All = NaN*ones(size(lon));
clme_AllF = clme_All;
clme_AllP = clme_All;
clme_AllD = clme_All;
pFracPD_grid = clme_All;
sFracPD_grid = clme_All;

for L=1:66
    lid = find(tlme==L);
    
    clme_All(lid) = plme_mcatch(L);
    clme_AllF(lid) = plme_Fmcatch(L);
    clme_AllP(lid) = plme_Pmcatch(L);
    clme_AllD(lid) = plme_Dmcatch(L);
    pFracPD_grid(lid) = pFracPD(L);
    sFracPD_grid(lid) = sFracPD(L);
end

diffP = pFracPD_grid - sFracPD_grid;

%% CATCH CORRELATIONS
figure(1)
subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10s(lme),l10p(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-1.75,1.7,['r = ' sprintf('%2.2f',rall)])
text(-1.75,1.4,['RMSE = ' sprintf('%2.2f',rmse)])
axis([-2 2 -2 2])
xlabel('SAUP')
ylabel('POEM')
title('D. All fishes')

subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10sF(lme),l10pF(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' sprintf('%2.2f',rF)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseF)])
axis([-6 2 -6 2])
xlabel('SAUP')
ylabel('POEM')
title('A. Forage Fishes')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10sP(lme),l10pP(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-5.5,1.5,['r = ' sprintf('%2.2f',rP)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseP)])
axis([-6 2 -6 2])
xlabel('SAUP')
ylabel('POEM')
title('B. Large Pelagics')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(l10sD(lme),l10pD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-1.75,1.7,['r = ' sprintf('%2.2f',rD)])
text(-1.75,1.4,['RMSE = ' sprintf('%2.2f',rmseD)])
axis([-2 2 -2 2])
xlabel('SAUP')
ylabel('POEM')
title('C. Demersals')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,cfile2,'_SAUP_comp_types_temp_Stock_LELC_ms.png'])

%% FRACTION LARGE PELAGIC Correlation
figure(2)
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(sFracPD(lme),pFracPD(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(sFracPD(lme),pFracPD(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center'); hold on;
end
text(0.05,0.95,['r = ' sprintf('%2.2f',rPD)])
text(0.05,0.9,['RMSE = ' sprintf('%2.2f',rmsePD)])
text(0.05,0.85,['Fmed = ' sprintf('%2.2f',FPD)])
axis([0 1 0 1])
xlabel('SAUP')
ylabel('POEM')
title('Fraction Large Pelagics')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,cfile2,'_SAUP_comp_fracP.png'])

%% Catch Map 
% all
% figure(3)
% clf
% subplot('Position',[0.5 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_All)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-3.5 1.5]);
% set(gcf,'renderer','painters')
% title('All Fishes')
% 
% % all F
% subplot('Position',[0 0.51 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_AllF)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-3.5 1.5]);
% %hcb = colorbar('h');
% colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
% set(gcf,'renderer','painters')
% title('mean log10 total annual catch (MT) Forage Fishes')
% 
% % all P
% subplot('Position',[0.5 0.51 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_AllP)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-3.5 1.5]);
% %hcb = colorbar('h');
% set(gcf,'renderer','painters')
% title('Large Pelagic Fishes')
% 
% % All D
% subplot('Position',[0 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_AllD)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-3.5 1.5]);
% %hcb = colorbar('h');
% set(gcf,'renderer','painters')
% title('Demersal Fishes')
% stamp(cfile)
% print('-dpng',[ppath 'Clim_fished_',harv,cfile2,'_LME_catch.png'])

%% subplot with SAU comparison map & corr
figure(4)
% POEM
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,pFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.01 0.5 0.5 0.05],'orientation','horizontal')                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology Fraction LP catch')

% Diff
subplot('Position',[0.0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('POEM - SAU difference')

% SAUP
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,sFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('SAU Fraction LP catch')

%Corr
subplot('Position',[0.55 0.075 0.4 0.4])
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(sFracPD(lme),pFracPD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.75,0.65,['r = ' sprintf('%2.2f',rPD)])
text(0.75,0.6,['RMSE = ' sprintf('%2.2f',rmsePD)])
text(0.75,0.55,['Fmed = ' sprintf('%2.2f',FPD)])
axis([0 1.05 0 1.05])
xlabel('SAUP')
ylabel('POEM')
title('Fraction Large Pelagics')
stamp(cfile)
print('-dpng',[ppath 'Clim_' harv cfile2 '_LME_fracPD_catch_SAUP_comp_subplot.png'])


