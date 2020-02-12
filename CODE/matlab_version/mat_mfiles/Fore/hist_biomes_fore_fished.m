% Calc biome biomass of FEISTY forecast
% MOM grid cells in 3 biomes
% LC: surf chl < 0.15 mg/m3
% ECCS: surf chl > 0.15 mg/m3; max MLD < 75 m
% ECSS: surf chl > 0.15 mg/m3; max MLD > 75 m
% Forecasted biomass on Hist biomes

clear all
close all

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/FEISTY/NC/Matlab_new_size/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([dpath 'Means_fore_',harv,'_' cfile '.mat']);
load([cpath 'COBALT_biomes_50yr_means_5100.mat'],'biome_hist');
load([gpath 'hindcast_gridspec.mat'],'geolat_t','geolon_t');
grid = csvread([gpath 'grid_csv.csv']);

%% Calc biome biomass avgs
gid = grid(:,1);

sf(gid) = sf_mean50(:,1);
sp(gid) = sp_mean50(:,1);
sd(gid) = sd_mean50(:,1);
mf(gid) = mf_mean50(:,1);
mp(gid) = mp_mean50(:,1);
md(gid) = md_mean50(:,1);
lp(gid) = lp_mean50(:,1);
ld(gid) = ld_mean50(:,1);
b(gid) = b_mean50(:,1);

FbiomeH_mbio50 = NaN*ones(3,9);
for L=1:3
    lid = find(biome_hist==L);
    FbiomeH_mbio50(L,1) = nanmean(sf(lid));
    FbiomeH_mbio50(L,2) = nanmean(sp(lid));
    FbiomeH_mbio50(L,3) = nanmean(sd(lid));
    FbiomeH_mbio50(L,4) = nanmean(mf(lid));
    FbiomeH_mbio50(L,5) = nanmean(mp(lid));
    FbiomeH_mbio50(L,6) = nanmean(md(lid));
    FbiomeH_mbio50(L,7) = nanmean(lp(lid));
    FbiomeH_mbio50(L,8) = nanmean(ld(lid));
    FbiomeH_mbio50(L,9) = nanmean(b(lid));
end

save([dpath 'Biomes_Fore_',harv,'_' cfile '.mat'],'FbiomeH_mbio50','-append');

%% Figures
FbiomeH_sf = NaN*ones(360,200);
FbiomeH_sp = FbiomeH_sf;
FbiomeH_sd = FbiomeH_sf;
FbiomeH_mf = FbiomeH_sf;
FbiomeH_mp = FbiomeH_sf;
FbiomeH_md = FbiomeH_sf;
FbiomeH_lp = FbiomeH_sf;
FbiomeH_ld = FbiomeH_sf;
FbiomeH_b = FbiomeH_sf;

for L=1:3
    lid = find(biome_hist==L);
    
    FbiomeH_sf(lid) = FbiomeH_mbio50(L,1);
    FbiomeH_sp(lid) = FbiomeH_mbio50(L,2);
    FbiomeH_sd(lid) = FbiomeH_mbio50(L,3);
    FbiomeH_mf(lid) = FbiomeH_mbio50(L,4);
    FbiomeH_mp(lid) = FbiomeH_mbio50(L,5);
    FbiomeH_md(lid) = FbiomeH_mbio50(L,6);
    FbiomeH_lp(lid) = FbiomeH_mbio50(L,7);
    FbiomeH_ld(lid) = FbiomeH_mbio50(L,8);
    FbiomeH_b(lid) = FbiomeH_mbio50(L,9);
end

FbiomeH_AllF = FbiomeH_sf+FbiomeH_mf;
FbiomeH_AllP = FbiomeH_sp+FbiomeH_mp+FbiomeH_lp;
FbiomeH_AllD = FbiomeH_sd+FbiomeH_md+FbiomeH_ld;
FbiomeH_All = FbiomeH_AllF + FbiomeH_AllP + FbiomeH_AllD;

%% Check map
% plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];
% ENTER -100 TO MAP ORIGIN LONG

%% ALL
figure(1)
subplot(2,2,4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(FbiomeH_All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
caxis([0 1])                    %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass All Fishes (g m^-^2)')
% stamp(cfile)
% print('-dpng',[ppath 'Fore_fished',harv,'_Histbiomes_All.png'])

% all F
subplot(2,2,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(FbiomeH_AllF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
caxis([-0.5 0.5])                 %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Forage Fishes (g m^-^2)')
% stamp(cfile)
% print('-dpng',[ppath 'Fore_fished',harv,'_Histbiomes_AllF.png'])

% all P
subplot(2,2,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(FbiomeH_AllP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
caxis([-0.5 0.5])                 %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Pelagic Fishes (g m^-^2)')
%stamp(cfile)
%print('-dpng',[ppath 'Fore_fished',harv,'_Histbiomes_AllP.png'])

% All D
subplot(2,2,3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(FbiomeH_AllD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
caxis([-0.5 0.5])               %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Demersal Fishes (g m^-^2)')
stamp(cfile)
%print('-dpng',[ppath 'Fore_fished',harv,'_Histbiomes_AllD.png'])
print('-dpng',[ppath 'Fore_fished',harv,'_Histbiomes_subplot.png'])

