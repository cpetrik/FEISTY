% Calc biome biomass of POEM
% MOM grid cells in 3 biomes
% LC: surf chl < 0.15 mg/m3
% ECCS: surf chl > 0.15 mg/m3; max MLD < 75 m
% ECSS: surf chl > 0.15 mg/m3; max MLD > 75 m

clear all
close all

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/FEISTY/NC/Matlab_new_size/';

cfile = 'Dc_enc70_cmax-metab20_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00400';
harv = '03';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([dpath 'Means_hist_fished',harv,'_' cfile '.mat']);
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

Hbiome_mbio50 = NaN*ones(3,9);
for L=1:3
    lid = find(biome_hist==L);
    Hbiome_mbio50(L,1) = nanmean(sf(lid));
    Hbiome_mbio50(L,2) = nanmean(sp(lid));
    Hbiome_mbio50(L,3) = nanmean(sd(lid));
    Hbiome_mbio50(L,4) = nanmean(mf(lid));
    Hbiome_mbio50(L,5) = nanmean(mp(lid));
    Hbiome_mbio50(L,6) = nanmean(md(lid));
    Hbiome_mbio50(L,7) = nanmean(lp(lid));
    Hbiome_mbio50(L,8) = nanmean(ld(lid));
    Hbiome_mbio50(L,9) = nanmean(b(lid));
end

save([dpath 'Biomes_hist_fished',harv,'_' cfile '.mat'],'Hbiome_mbio50');


%% Figures
Hbiome_sf = NaN*ones(360,200);
Hbiome_sp = Hbiome_sf;
Hbiome_sd = Hbiome_sf;
Hbiome_mf = Hbiome_sf;
Hbiome_mp = Hbiome_sf;
Hbiome_md = Hbiome_sf;
Hbiome_lp = Hbiome_sf;
Hbiome_ld = Hbiome_sf;
Hbiome_b = Hbiome_sf;

for L=1:3
    lid = find(biome_hist==L);
    
    Hbiome_sf(lid) = Hbiome_mbio50(L,1);
    Hbiome_sp(lid) = Hbiome_mbio50(L,2);
    Hbiome_sd(lid) = Hbiome_mbio50(L,3);
    Hbiome_mf(lid) = Hbiome_mbio50(L,4);
    Hbiome_mp(lid) = Hbiome_mbio50(L,5);
    Hbiome_md(lid) = Hbiome_mbio50(L,6);
    Hbiome_lp(lid) = Hbiome_mbio50(L,7);
    Hbiome_ld(lid) = Hbiome_mbio50(L,8);
    Hbiome_b(lid) = Hbiome_mbio50(L,9);
end

Hbiome_AllF = Hbiome_sf+Hbiome_mf;
Hbiome_AllP = Hbiome_sp+Hbiome_mp+Hbiome_lp;
Hbiome_AllD = Hbiome_sd+Hbiome_md+Hbiome_ld;
Hbiome_All = Hbiome_AllF + Hbiome_AllP + Hbiome_AllD;

%% Check map
% plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
% ENTER -100 TO MAP ORIGIN LONG

% Pac only
lonlim=[-255 -60];
latlim=[0 80];

%% ALL
figure(41)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(Hbiome_All)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.3 0.9]);
hcb = colorbar('h');
ylim(hcb,[0.3 0.9])                    %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass All Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished',harv,'_Pac_biomes_All.png'])

% all F
figure(42)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(Hbiome_AllF)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.8 0.8]);
hcb = colorbar('h');
ylim(hcb,[-0.8 0.8])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Forage Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished',harv,'_Pac_biomes_AllF.png'])

% all P
figure(43)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(Hbiome_AllP)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.8 0.8]);
hcb = colorbar('h');
ylim(hcb,[-0.8 0.8])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Pelagic Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished',harv,'_Pac_biomes_AllP.png'])

% All D
figure(44)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(Hbiome_AllD)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.8 0.8]);
hcb = colorbar('h');
ylim(hcb,[-0.8 0.8])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Demersal Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished',harv,'_Pac_biomes_AllD.png'])
