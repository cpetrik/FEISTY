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
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

cfile = 'Dc_enc70_cmax-metab20_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00400';
harv = '03';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([dpath 'Means_fore_pristine_' cfile '.mat']);
load([cpath 'COBALT_biomes_last50yr.mat'],'biome_fore');
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

biome_mbio50 = NaN*ones(3,9);
for L=1:3
    lid = find(biome_fore==L);
    biome_mbio50(L,1) = nanmean(sf(lid));
    biome_mbio50(L,2) = nanmean(sp(lid));
    biome_mbio50(L,3) = nanmean(sd(lid));
    biome_mbio50(L,4) = nanmean(mf(lid));
    biome_mbio50(L,5) = nanmean(mp(lid));
    biome_mbio50(L,6) = nanmean(md(lid));
    biome_mbio50(L,7) = nanmean(lp(lid));
    biome_mbio50(L,8) = nanmean(ld(lid));
    biome_mbio50(L,9) = nanmean(b(lid));
end

save([dpath 'Biomes_fore_pristine_' cfile '.mat'],'biome_mbio50');


%% Figures
biome_sf = NaN*ones(360,200);
biome_sp = biome_sf;
biome_sd = biome_sf;
biome_mf = biome_sf;
biome_mp = biome_sf;
biome_md = biome_sf;
biome_lp = biome_sf;
biome_ld = biome_sf;
biome_b = biome_sf;

for L=1:3
    lid = find(biome_fore==L);
    
    biome_sf(lid) = biome_mbio50(L,1);
    biome_sp(lid) = biome_mbio50(L,2);
    biome_sd(lid) = biome_mbio50(L,3);
    biome_mf(lid) = biome_mbio50(L,4);
    biome_mp(lid) = biome_mbio50(L,5);
    biome_md(lid) = biome_mbio50(L,6);
    biome_lp(lid) = biome_mbio50(L,7);
    biome_ld(lid) = biome_mbio50(L,8);
    biome_b(lid) = biome_mbio50(L,9);
end

biome_All = biome_sf+biome_sp+biome_sd+biome_mf+biome_mp+biome_md+biome_lp+biome_ld;
biome_AllF = biome_sf+biome_mf;
biome_AllP = biome_sp+biome_mp+biome_lp;
biome_AllD = biome_sd+biome_md+biome_ld;

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
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(biome_All)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.3 0.9]);
hcb = colorbar('h');
ylim(hcb,[0.3 0.9])                    %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass All Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Fore_pristine_biomes_All.png'])

%% all F
figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(biome_AllF)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.8 0.8]);
hcb = colorbar('h');
ylim(hcb,[-0.8 0.8])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Forage Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Fore_pristine_biomes_AllF.png'])

% all P
figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(biome_AllP)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.8 0.8]);
hcb = colorbar('h');
ylim(hcb,[-0.8 0.8])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Pelagic Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Fore_pristine_biomes_AllP.png'])

% All D
figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(biome_AllD)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.8 0.8]);
hcb = colorbar('h');
ylim(hcb,[-0.8 0.8])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Demersal Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Fore_pristine_biomes_AllD.png'])

%% Pac only
lonlim=[-255 -60];
latlim=[0 80];

% ALL
figure(41)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(biome_All)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.3 0.9]);
hcb = colorbar('h');
ylim(hcb,[0.3 0.9])                    %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass All Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Fore_pristine_Pac_biomes_All.png'])

% all F
figure(42)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(biome_AllF)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.8 0.8]);
hcb = colorbar('h');
ylim(hcb,[-0.8 0.8])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Forage Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Fore_pristine_Pac_biomes_AllF.png'])

% all P
figure(43)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(biome_AllP)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.8 0.8]);
hcb = colorbar('h');
ylim(hcb,[-0.8 0.8])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Pelagic Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Fore_pristine_Pac_biomes_AllP.png'])

% All D
figure(44)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(biome_AllD)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.8 0.8]);
hcb = colorbar('h');
ylim(hcb,[-0.8 0.8])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Demersal Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Fore_pristine_Pac_biomes_AllD.png'])
