% P:D ratio by LME
% Climatology
% 150 years
% Saved as mat files
% Compare to Daniel's model results & Reg

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'LME_clim_temp_zoop_det.mat']);

tlme = lme_mask_onedeg;
AREA_OCN = max(area,1);

% "lfile" never changes, has lme areas
lfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
lpath = ['/Volumes/GFDL/NC/Matlab_new_size/' lfile '/'];
load([lpath 'LME_clim_fished03_' lfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

%% DvD on grid
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/DanielVD_PelDem/Colleen_modeledfish_LME.mat')
FracLP_grid = NaN*ones(180,360);
for L=1:63
    lid = find(tlme==L);
    FracLP_grid(lid) = FracLP(L);
end

%% Reg
%use weighted catches
load([spath 'Reg_LME_Catch_annual.mat']);

% sum of functional types in each
Plme_land_all = nansum(Plme_wland,3);
Dlme_land_all = nansum(Dlme_wland,3);
Plme_catch_all = nansum(Plme_wall,3);
Dlme_catch_all = nansum(Dlme_wall,3);

% 1970-2014 Reg average
%Mean of that time period
Plme_mland10 = nanmean(Plme_land_all)';
Dlme_mland10 = nanmean(Dlme_land_all)';
Plme_mcatch10 = nanmean(Plme_catch_all)';
Dlme_mcatch10 = nanmean(Dlme_catch_all)';

% MT/km2
Plme_mland10 = Plme_mland10 ./ lme_area_km2;
Dlme_mland10 = Dlme_mland10 ./ lme_area_km2;
Plme_mcatch10 = Plme_mcatch10 ./ lme_area_km2;
Dlme_mcatch10 = Dlme_mcatch10 ./ lme_area_km2;

lFracPD = Plme_mland10 ./ (Plme_mland10 + Dlme_mland10);
cFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

%on grid
tlme = lme_mask_onedeg;
lFracPD_grid = NaN*ones(180,360);
cFracPD_grid = NaN*ones(180,360);
for L=1:66
    lid = find(tlme==L);
    lFracPD_grid(lid) = lFracPD(L);
    cFracPD_grid(lid) = cFracPD(L);
end

%% Comparison stats
keep=[1:61,63];
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
notLELC = notLELC(notLELC<=63);

diffL = FracLP_grid - lFracPD_grid;
diffC = FracLP_grid - cFracPD_grid;

%r
rL=corr(FracLP(keep),lFracPD(keep));
rC=corr(FracLP(keep),cFracPD(keep));

%root mean square error
o=FracLP(keep);
p=lFracPD(keep);
n = length(o);
num=nansum((p-o).^2);
rmseL = sqrt(num/n);

o=FracLP(keep);
p=cFracPD(keep);
n = length(o);
num=nansum((p-o).^2);
rmseC = sqrt(num/n);

%Fmed
FL=10^(median(FracLP(keep)-lFracPD(keep)));
FC=10^(median(FracLP(keep)-cFracPD(keep)));

% Table
fish_stat(1,1) = rL;
fish_stat(2,1) = rmseL;
fish_stat(3,1) = FL;
fish_stat(1,2) = rC;
fish_stat(2,2) = rmseC;
fish_stat(3,2) = FC;

Fstat = array2table(fish_stat,'RowNames',{'r','RMSE','Fmed'},...
    'VariableNames',{'Landings','Catch'});
writetable(Fstat,[spath 'LME_DvDyrs_vs_Reg_stats.csv'],'Delimiter',',','WriteRowNames',true)
save([spath 'LME_DvDyrs_vs_Reg_stats.mat'],'fish_stat')

%% Plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
% ENTER -100 TO MAP ORIGIN LONG

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

revamp=colormap(cmocean('amp'));
revamp=flipud(revamp);
% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

x=0:0.1:1;

%% Correlation
figure(1)
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(lFracPD(lme),FracLP(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(lFracPD(lme),FracLP(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center'); hold on;
end
text(0.05,0.95,['r = ' sprintf('%2.2f',rL)])
text(0.05,0.9,['RMSE = ' sprintf('%2.2f',rmseL)])
text(0.05,0.85,['Fmed = ' sprintf('%2.2f',FL)])
axis([-0.05 1.05 -0.05 1.05])
xlabel('Reg landings')
ylabel('DvD')
title('Fraction Large Pelagics')
stamp('')
print('-dpng',[spath 'LME_DvDyrs_vs_Reg_landings_comp_fracP.png'])

figure(2)
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(cFracPD(lme),FracLP(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(cFracPD(lme),FracLP(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center'); hold on;
end
text(0.05,0.95,['r = ' sprintf('%2.2f',rC)])
text(0.05,0.9,['RMSE = ' sprintf('%2.2f',rmseC)])
text(0.05,0.85,['Fmed = ' sprintf('%2.2f',FC)])
axis([-0.05 1.05 -0.05 1.05])
xlabel('Reg catch')
ylabel('DvD')
title('Fraction Large Pelagics')
stamp('')
print('-dpng',[spath 'LME_DvDyrs_vs_Reg_catch_comp_fracP.png'])


%% Subplot with maps and corr
figure(3)
% DvD
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracLP_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.01 0.5 0.5 0.05],'orientation','horizontal')                   %Set color axis if needed
set(gcf,'renderer','painters')
title('DvD Fraction LP catch')

% Diff
subplot('Position',[0.0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('DvD - Reg difference')

%Reg
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Reg Fraction LP landings')

%Corr
subplot('Position',[0.55 0.075 0.4 0.4])
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(lFracPD(lme),FracLP(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.05,0.95,['r = ' sprintf('%2.2f',rL)])
text(0.05,0.9,['RMSE = ' sprintf('%2.2f',rmseL)])
text(0.05,0.85,['Fmed = ' sprintf('%2.2f',FL)])
axis([0 1.05 0 1.05])
xlabel('Reg landings')
ylabel('DvD')
title('Fraction Large Pelagics')
stamp('')
print('-dpng',[spath 'LME_fracPD_DvDyrs_vs_Reg_landings_comp_subplot.png'])


figure(4)
% DvD
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracLP_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.01 0.5 0.5 0.05],'orientation','horizontal')                   %Set color axis if needed
set(gcf,'renderer','painters')
title('DvD Fraction LP catch')

% Diff
subplot('Position',[0.0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffC)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('DvD - Reg difference')

%Reg
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,cFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Reg Fraction LP catch')

%Corr
subplot('Position',[0.55 0.075 0.4 0.4])
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(cFracPD(lme),FracLP(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.05,0.95,['r = ' sprintf('%2.2f',rL)])
text(0.05,0.9,['RMSE = ' sprintf('%2.2f',rmseL)])
text(0.05,0.85,['Fmed = ' sprintf('%2.2f',FL)])
axis([0 1.05 0 1.05])
xlabel('Reg catch')
ylabel('DvD')
title('Fraction Large Pelagics')
stamp('')
print('-dpng',[spath 'LME_fracPD_DvDyrs_vs_Reg_catch_comp_subplot.png'])



