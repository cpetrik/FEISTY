% P:D ratio by LME 
% Climatology
% 150 years
% Saved as mat files
% Compare to Daniel's model results &SAUP

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
load([cdir 'temp_100_1deg_ESM26_5yr_clim_191_195.mat'])
load([cdir 'btm_temp_1deg_ESM26_5yr_clim_191_195.mat'])
load([cpath 'LME_clim_temp_zoop_det.mat']);

ptemp_mean_clim=squeeze(nanmean(temp_100,1));
btemp_mean_clim=squeeze(nanmean(btm_temp,1));
tlme = lme_mask_onedeg;
AREA_OCN = max(area,1);

% "lfile" never changes, has lme areas
lfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
lpath = ['/Volumes/GFDL/NC/Matlab_new_size/' lfile '/'];
load([lpath 'LME_clim_fished03_' lfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

%% DvD on grid
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/DanielVD_PelDem/Colleen_modeledfish_LME.mat')
dlme_Pfrac = NaN*ones(180,360);
for L=1:63
    lid = find(tlme==L);
    dlme_Pfrac(lid) = FracLP(L);
end

%% SAUP
%use weighted catches
load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
    'Dlme_wcatch','Plme_wcatch');

Plme_catch_all = nansum(Plme_wcatch,3);
Dlme_catch_all = nansum(Dlme_wcatch,3);
lme_catch_all = nansum(lme_catch,3);

%1950-2006 SAUP average
id = find(yr>1950 & yr<=2006);

slme_mcatch = nanmean(lme_catch_all(id,:));
slme_mcatch = slme_mcatch';
Pslme_mcatch = nanmean(Plme_catch_all(id,:));
Pslme_mcatch = Pslme_mcatch';
Dslme_mcatch = nanmean(Dlme_catch_all(id,:));
Dslme_mcatch = Dslme_mcatch';

slme_mcatch10 = NaN*ones(size(slme_mcatch));
Plme_mcatch10 = NaN*ones(size(slme_mcatch));
Dlme_mcatch10 = NaN*ones(size(slme_mcatch));

%Top 10 yrs by LME SAUP 
for i=1:66
    [sort_lme_catch,ix] = sort(lme_catch_all(:,i),'descend');
    sort_Plme_catch = Plme_catch_all(ix,i);
    sort_Dlme_catch = Dlme_catch_all(ix,i);
    slme_mcatch10(i) = nanmean(sort_lme_catch(1:10));
    Plme_mcatch10(i) = nanmean(sort_Plme_catch(1:10));
    Dlme_mcatch10(i) = nanmean(sort_Dlme_catch(1:10));
end

% MT/km2
Plme_mcatch10 = Plme_mcatch10 ./ lme_area_km2;
Dlme_mcatch10 = Dlme_mcatch10 ./ lme_area_km2;

sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

%on grid
tlme = lme_mask_onedeg;
sFracPD_grid = NaN*ones(180,360);
for L=1:66
    lid = find(tlme==L);
    sFracPD_grid(lid) = sFracPD(L);
end

%% Comparison stats
keep=[1:61,63];
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
notLELC = notLELC(notLELC<=63);

diffD = dlme_Pfrac - sFracPD_grid;

%r
rall=corr(FracLP(keep),sFracPD(keep));
r2=corr(FracLP(notLELC),sFracPD(notLELC));

%root mean square error
o=FracLP(keep);
p=sFracPD(keep);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=FracLP(notLELC);
p=sFracPD(notLELC);
n = length(o);
num=nansum((p-o).^2);
rmse2 = sqrt(num/n);

%Fmed
Fall=10^(median(FracLP(keep)-sFracPD(keep)));
Fall2=10^(median(FracLP(notLELC)-sFracPD(notLELC)));

% Table
fish_stat(1,1) = rall;
fish_stat(2,1) = rmse;
fish_stat(3,1) = Fall;
fish_stat(1,2) = r2;
fish_stat(2,2) = rmse2;
fish_stat(3,2) = Fall2;

Fstat = array2table(fish_stat,'RowNames',{'r','RMSE','Fmed'},...
    'VariableNames',{'AllLMEs','noLELC'});
writetable(Fstat,[spath 'LME_DvD_vs_SAU_stats.csv'],'Delimiter',',','WriteRowNames',true)
save([spath 'LME_DvD_vs_SAU_stats.mat'],'fish_stat')

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
    plot(sFracPD(lme),FracLP(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(sFracPD(lme),FracLP(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center'); hold on;
end
text(0.05,0.95,['r = ' sprintf('%2.2f',rall)])
text(0.05,0.9,['RMSE = ' sprintf('%2.2f',rmse)])
text(0.05,0.85,['Fmed = ' sprintf('%2.2f',Fall)])
axis([-0.05 1.05 -0.05 1.05])
xlabel('SAUP')
ylabel('DvD')
title('Fraction Large Pelagics')
stamp('')
print('-dpng',[spath 'LME_DvD_vs_SAUP_comp_fracP.png'])


%% Subplot with maps and corr
figure(2)
% DvD
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,dlme_Pfrac)
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
surfm(geolat_t,geolon_t,diffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('DvD - SAU difference')

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
    plot(sFracPD(lme),FracLP(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.05,0.95,['r = ' sprintf('%2.2f',rall)])
text(0.05,0.9,['RMSE = ' sprintf('%2.2f',rmse)])
text(0.05,0.85,['Fmed = ' sprintf('%2.2f',Fall)])
axis([0 1.05 0 1.05])
xlabel('SAUP')
ylabel('DvD')
title('Fraction Large Pelagics')
stamp('')
print('-dpng',[spath 'LME_fracPD_DvD_vs_SAUP_comp_subplot.png'])



