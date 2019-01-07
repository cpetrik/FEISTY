%POEM TEs vs. Mauread ECIs by LME

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'LME_clim_temp_zoop_det.mat']);

load([spath 'Maureaud_etal_2017_s002_ECI.mat']);

% POEM file info
frate = 0.3;
tfish = num2str(100+int64(10*frate));

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BEdecrT03_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'TEeff_Climatol_All_fish03_' cfile '.mat']);

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')

%% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

x=-6:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% ECI for clim years (1991-1995?)
mECI = mean(ECI(:,2:6),2);
keep = ECI(:,1);

%% POEM LME TEeffs
%     lme_te(L,2) = nanmean(TEeff_L(lid));
%     lme_te(L,4) = nanmean(TEeff_HTLd(lid));
%     lme_te(L,6) = nanmean(TEeff_LTLd(lid));
pECI = lme_te(keep,4);

%% Drop same as Mauread -------------------------
% table
tab(:,1)=keep;
tab(:,2)=mECI;
tab(:,3)=pECI;

T1 = array2table(tab,'VariableNames',{'LME','Mauread','POEM'});

writetable(T1,[spath 'LME_TEeff_Mauread_comp_' cfile '.csv'],'Delimiter',',')
save([dpath 'LME_TEeff_Mauread_comp_' cfile '.mat'],'tab','keep',...
    'mECI','pECI')

%% Stats
ma = (mECI); 
po = (pECI); 
Lma = log10(mECI); %log10(mECI)
Lpo = log10(pECI);

%r
rall=corr(ma,po);
rL=corr(Lma,Lpo);

%root mean square error
o=ma;
p=po;
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=Lma;
p=Lpo;
n = length(o);
num=nansum((p-o).^2);
rmseL = sqrt(num/n);

%Fmed
Fall=10^(median(ma-po));
FL=10^(median(Lma-Lpo));

% Table
fish_stat(1,1) = rall;
fish_stat(1,2) = rmse;
fish_stat(1,3) = Fall;
fish_stat(2,1) = rL;
fish_stat(2,2) = rmseL;
fish_stat(2,3) = FL;

Fstat = array2table(fish_stat,'VariableNames',{'r','RMSE','Fmed'},...
    'RowNames',{'raw','log10'});
writetable(Fstat,[dpath 'LME_TEeff_Mauread_stats_' cfile '.csv'],'Delimiter',',',...
    'WriteRowNames',true)
save([dpath 'LME_TEeff_Mauread_stats_' cfile '.mat'],'fish_stat')

%% Plots
figure(1)
plot(x,x,'--k'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(mECI(i),pECI(i),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.035,0.03,['r = ' sprintf('%2.2f',rall)])
text(0.035,0.0275,['RMSE = ' sprintf('%2.2f',rmse)])
text(0.035,0.025,['Fmed = ' sprintf('%2.2f',Fall)])
axis([0 0.045 0 0.045])
xlabel('Mauread ECI')
ylabel('POEM TEeff HTL')
title('Transfer efficiency by LME')
print('-dpng',[ppath 'LME_TEeff_Mauread_corr_raw.png']);

figure(2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':k'); hold on;
plot(x,x2l,':k'); hold on;
% plot(x,x5h,'--r'); hold on;
% plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(log10(mECI(i)),log10(pECI(i)),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-1.5,-2.1,['r = ' sprintf('%2.2f',rL)])
text(-1.5,-2.2,['RMSE = ' sprintf('%2.2f',rmseL)])
text(-1.5,-2.3,['Fmed = ' sprintf('%2.2f',FL)])
axis([-3.5 -1 -3.5 -1])
xlabel('Mauread ECI')
ylabel('POEM TEeff HTL')
title('log10 Transfer efficiency by LME')
print('-dpng',[ppath 'LME_TEeff_Mauread_corr_log10.png']);

figure(3)
subplot(1,2,1)
plot(x,x,'--k'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(mECI(i),pECI(i),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.025,0.01,['r = ' sprintf('%2.2f',rall)])
text(0.025,0.0075,['RMSE = ' sprintf('%2.2f',rmse)])
text(0.025,0.005,['Fmed = ' sprintf('%2.2f',Fall)])
axis([0 0.045 0 0.045])
axis equal
xlabel('Mauread ECI')
ylabel('POEM TEeff HTL')
title('Transfer efficiency by LME')

subplot(1,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':k'); hold on;
plot(x,x2l,':k'); hold on;
% plot(x,x5h,'--r'); hold on;
% plot(x,x5l,'--r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(log10(mECI(i)),log10(pECI(i)),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-1.95,-3.05,['r = ' sprintf('%2.2f',rL)])
text(-1.95,-3.2,['RMSE = ' sprintf('%2.2f',rmseL)])
text(-1.95,-3.35,['Fmed = ' sprintf('%2.2f',FL)])
axis([-3.5 -1 -3.5 -1])
axis equal
xlabel('Mauread ECI')
ylabel('POEM TEeff HTL')
title('log10 Transfer efficiency by LME')
print('-dpng',[ppath 'LME_TEeff_Mauread_corr_subplot.png']);


%% Catch map
%on grid
tlme = lme_mask_onedeg;
mECI_grid = NaN*ones(180,360);
pECI_grid = NaN*ones(180,360);
for L=1:length(keep)
    lme=keep(L);
    lid = find(tlme==lme);
    mECI_grid(lid) = mECI(L);
    pECI_grid(lid) = pECI(L);
end


diffP = pECI_grid - mECI_grid;
%%
figure(4)
% POEM
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,pECI_grid)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.005 0.035]);
colorbar('Position',[0.025 0.555 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('Climatology TEeff HTL')

% Diff
subplot('Position',[0.25 0.025 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.03 0.03]);
colorbar('Position',[0.275 0.05 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('POEM - Mauread difference')

% Mauread
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,mECI_grid)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.005 0.035]);
colorbar('Position',[0.525 0.555 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('Mauread mean ECI 1991-1995')
stamp(cfile)
print('-dpng',[ppath 'Clim_' harv '_LME_TEeff_Mauread_comp.png'])

%% log10
diffL = log10(pECI_grid) - log10(mECI_grid);
figure(4)
% POEM
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(pECI_grid))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 -1.5]);
colorbar('Position',[0.025 0.555 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('Climatology log10 TEeff HTL')

% Diff
subplot('Position',[0.25 0.025 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.275 0.05 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('POEM - Mauread log10 difference')

% Mauread
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(mECI_grid))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 -1.5]);
colorbar('Position',[0.525 0.555 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('Mauread log10 mean ECI 1991-1995')
stamp(cfile)
print('-dpng',[ppath 'Clim_' harv '_LME_TEeff_Mauread_comp_log10.png'])

%% Subplot with maps and corr
figure(6)
% POEM
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(pECI_grid))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 -1]);
colorbar('Position',[0.025 0.555 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('Climatology log10 TEeff HTL')
text(-2.75,1.25,'A')

% Diff
subplot('Position',[0.0 0.025 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.025 0.05 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('POEM - Mauread log10 difference')
text(-2.75,1.25,'C')

% Mauread
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(mECI_grid))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 -1]);
colorbar('Position',[0.525 0.555 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('Mauread log10 mean ECI 1991-1995')
text(-2.75,1.25,'B')

%Corr
subplot('Position',[0.575 0.075 0.375 0.375])
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(log10(mECI(i)),log10(pECI(i)),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-1.6,-2.2,['r = ' sprintf('%2.2f',rL)])
text(-1.6,-2.4,['RMSE = ' sprintf('%2.2f',rmseL)])
%text(-1.5,-3.3,['Fmed = ' sprintf('%2.2f',FL)])
axis([-3 -1 -3 -1])
xlabel('Mauread ECI')
ylabel('POEM TEeff HTL')
title('log10 Transfer efficiency by LME')
text(-2.9,-1.15,'D')
%stamp(cfile)
print('-dpng',[ppath 'Clim_' harv '_LME_TEeff_Mauread_comp_subplot.png'])

%% TEST LOCATIONS ============================================
locs=[8;19;22;21;20;1;51;49;10;6;13];
name = {'SS','GS','NS','NwS','BS','EBS','OyaCur','KurCur','Haw','SEUS',...
    'Humb'};

Tab=table(locs,name',mECI(locs),pECI(locs),...
    'VariableNames',{'LME','Name','Mauread','POEM'});
writetable(Tab,[spath 'LME_locs_TEeff_Mauread_',harv,'_' cfile '.csv'],'Delimiter',',');
save([dpath 'LME_locs_TEeff_Mauread_',harv,'_' cfile '.mat'],'Tab');



