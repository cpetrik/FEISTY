%FEISTY  TEs vs. Maureaud ECIs by LME

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

% FEISTY  file info
frate = 0.3;
tfish = num2str(100+int64(10*frate));

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
BE = 0.075;
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'TEeff_Climatol_All_fish03_' cfile '.mat']);

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')
warning off
cmYOR=cbrewer('seq','YlOrRd',66);
cmRP=cbrewer('seq','RdPu',28);
cmPR=cbrewer('seq','PuRd',28);


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

x=-6:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% ECI for clim years (1991-1995?)
mECI = mean(ECI(:,2:6),2);
keep = ECI(:,1);

%% FEISTY  LME TEeffs
%     lme_te(L,2) = nanmean(TEeff_L(lid));
%     lme_te(L,4) = nanmean(TEeff_HTLd(lid));
%     lme_te(L,6) = nanmean(TEeff_LTLd(lid));
pECI = lme_te(keep,4);

%% Drop same as Maureaud -------------------------
% table
tab(:,1)=keep;
tab(:,2)=mECI;
tab(:,3)=pECI;

T1 = array2table(tab,'VariableNames',{'LME','Maureaud','POEM'});

writetable(T1,[spath 'LME_TEeff_Mauread_comp_' cfile '.csv'],'Delimiter',',')
save([dpath 'LME_TEeff_Mauread_comp_' cfile '.mat'],'tab','keep',...
    'mECI','pECI')

%% Stats
ma = (mECI); 
po = (pECI); 
Lma = log10(mECI); %log10(mECI)
Lpo = log10(pECI);

%r
[rall,pall]=corr(ma,po);
[rL,pL]=corr(Lma,Lpo);

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

% log10
diffL = log10(pECI_grid) - log10(mECI_grid);

%% Subplot with maps and corr
figure(1)
% POEM
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(pECI_grid))
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2.8 -1.4]);
colorbar('Position',[0.025 0.555 0.45 0.035],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('Climatology log_1_0 TEeff HTL')
text(-2.75,1.25,'A')

% Maureaud
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(mECI_grid))
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2.8 -1.4]);
colorbar('Position',[0.525 0.555 0.45 0.035],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('Maureaud log_1_0 mean ECI 1991-1995')
text(-2.75,1.25,'B')

% Diff
subplot('Position',[0.0 0.03 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.025 0.05 0.45 0.035],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('FEISTY  - Maureaud log_1_0 difference')
text(-2.75,1.25,'C')

%Corr
subplot('Position',[0.58 0.115 0.34 0.34])
plot(x,x,'--k');hold on;
scatter(log10(mECI),log10(pECI),20,lme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.935 0.11 0.025 0.34])
text(-1.65,-2.2,['r = ' sprintf('%2.2f',rL)])
text(-1.65,-2.4,['(p = ' sprintf('%2.2f',pL) ')'])
text(-1.65,-2.6,['RMSE = ' sprintf('%2.2f',rmseL)])
axis([-3 -1 -3 -1])
xlabel('Maureaud ECI')
ylabel('FEISTY  TEeff HTL')
title('log_1_0 Transfer efficiency by LME')
text(-2.9,-1.15,'D')
%stamp(cfile)
print('-dpng',[ppath 'FigS6_Clim_' harv '_LME_TEeff_Maureaud_comp_subplot_pval.png'])



