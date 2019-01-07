% Calc LME biomass of POEM
% Forecast time period 2006-2100
% 95 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

cfile = 'Dc_enc70_cmax-metab20_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00400';
harv = '03';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'Means_fore_pristine_' cfile '.mat']);

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t','AREA_OCN');
grid = csvread([cpath 'grid_csv.csv']);
load([cpath 'lme_mask_esm2m.mat']);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

%% Plots in space
[ni,nj]=size(geolon_t);

Zmf=NaN*ones(ni,nj);
Bsf = Zmf;
Bsp = Zmf;
Bsd = Zmf;
Bmf = Zmf;
Bmp = Zmf;
Bmd = Zmf;
Blp = Zmf;
Bld = Zmf;
Bb = Zmf;

Bsf(grid(:,1))=sf_mean50;
Bsp(grid(:,1))=sp_mean50;
Bsd(grid(:,1))=sd_mean50;
Bmf(grid(:,1))=mf_mean50;
Bmp(grid(:,1))=mp_mean50;
Bmd(grid(:,1))=md_mean50;
Blp(grid(:,1))=lp_mean50;
Bld(grid(:,1))=ld_mean50;
Bb(grid(:,1))=b_mean50;

% g/m2 --> total g
Asf_mean = Bsf .* AREA_OCN * 365;
Asp_mean = Bsp .* AREA_OCN * 365;
Asd_mean = Bsd .* AREA_OCN * 365;
Amf_mean = Bmf .* AREA_OCN * 365;
Amp_mean = Bmp .* AREA_OCN * 365;
Amd_mean = Bmd .* AREA_OCN * 365;
Alp_mean = Blp .* AREA_OCN * 365;
Ald_mean = Bld .* AREA_OCN * 365;
Ab_mean = Bb .* AREA_OCN * 365;

%% Calc LMEs
lat = geolat_t;
lon = geolon_t;

tlme = lme_mask_esm2m';

lme_mbio = NaN*ones(66,9);

for L=1:66
    lid = find(tlme==L);
    %mean biomass
    lme_mbio(L,1) = nanmean(Asf_mean(lid));
    lme_mbio(L,2) = nanmean(Asp_mean(lid));
    lme_mbio(L,3) = nanmean(Asd_mean(lid));
    lme_mbio(L,4) = nanmean(Amf_mean(lid));
    lme_mbio(L,5) = nanmean(Amp_mean(lid));
    lme_mbio(L,6) = nanmean(Amd_mean(lid));
    lme_mbio(L,7) = nanmean(Alp_mean(lid));
    lme_mbio(L,8) = nanmean(Ald_mean(lid));
    lme_mbio(L,9) = nanmean(Ab_mean(lid));
end

%%
save([dpath 'LME_fore_pristine_' cfile '.mat'],...
    'lme_mbio');

%% Figures

lme_sf = NaN*ones(360,200);
lme_sp = lme_sf;
lme_sd = lme_sf;
lme_mf = lme_sf;
lme_mp = lme_sf;
lme_md = lme_sf;
lme_lp = lme_sf;
lme_ld = lme_sf;
lme_b = lme_sf;

for L=1:66
    lid = find(tlme==L);
    
    lme_sf(lid) = lme_mbio(L,1);
    lme_sp(lid) = lme_mbio(L,2);
    lme_sd(lid) = lme_mbio(L,3);
    lme_mf(lid) = lme_mbio(L,4);
    lme_mp(lid) = lme_mbio(L,5);
    lme_md(lid) = lme_mbio(L,6);
    lme_lp(lid) = lme_mbio(L,7);
    lme_ld(lid) = lme_mbio(L,8);
    lme_b(lid) = lme_mbio(L,9);
end

lme_All = lme_sf+lme_sp+lme_sd+lme_mf+lme_mp+lme_md+lme_lp+lme_ld;
lme_AllF = lme_sf+lme_mf;
lme_AllP = lme_sp+lme_mp+lme_lp;
lme_AllD = lme_sd+lme_md+lme_ld;

% plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
% ENTER -100 TO MAP ORIGIN LONG

land=-999*ones(ni,nj);
land(grid(:,1))=NaN*ones(size(mf_mean50));

