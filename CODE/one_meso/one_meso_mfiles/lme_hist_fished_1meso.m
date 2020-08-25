% Calc LME biomass of FEISTY
% Hindcast time period 1951-2000 (for comparing with RCP 8.5)
% 50 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/FEISTY/NC/Matlab_new_size/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t','AREA_OCN');
grid = csvread([cpath 'grid_csv.csv']);
load([cpath 'lme_mask_esm2m.mat']);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

%% FEISTY
%cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100_noHPloss';
cfile = 'Dc_Lam700_enc70-b250_m400-b175-k086_c20-b250_D080_A050_nmort1_BE08_CC80_RE00100';
harv = 'All_fish03';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];
% dpath = [dp cfile '_noHPloss/'];
% ppath = [pp cfile '_noHPloss/'];

load([dpath 'Means_Historic_1meso_',harv,'_' cfile '.mat']);

%% Plots in space
[ni,nj]=size(geolon_t);

Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);

Zmf(grid(:,1))=mf_my50;
Zmp(grid(:,1))=mp_my50;
Zmd(grid(:,1))=md_my50;
Zlp(grid(:,1))=lp_my50;
Zld(grid(:,1))=ld_my50;

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
Amf_mcatch = Zmf .* AREA_OCN * 365; %mean fish catch per yr
Amp_mcatch = Zmp .* AREA_OCN * 365;
Amd_mcatch = Zmd .* AREA_OCN * 365;
Alp_mcatch = Zlp .* AREA_OCN * 365;
Ald_mcatch= Zld .* AREA_OCN * 365;

Asf_mean = Bsf .* AREA_OCN;
Asp_mean = Bsp .* AREA_OCN;
Asd_mean = Bsd .* AREA_OCN;
Amf_mean = Bmf .* AREA_OCN;
Amp_mean = Bmp .* AREA_OCN;
Amd_mean = Bmd .* AREA_OCN;
Alp_mean = Blp .* AREA_OCN;
Ald_mean = Bld .* AREA_OCN;
Ab_mean = Bb .* AREA_OCN;

%% Calc LMEs
lat = geolat_t;
lon = geolon_t;

tlme = lme_mask_esm2m';

lme_mcatch = NaN*ones(66,5);
lme_mbio = NaN*ones(66,9);
lme_sbio = NaN*ones(66,9);
lme_area = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    %total catch g
    lme_mcatch(L,1) = nansum(Amf_mcatch(lid));
    lme_mcatch(L,2) = nansum(Amp_mcatch(lid));
    lme_mcatch(L,3) = nansum(Amd_mcatch(lid));
    lme_mcatch(L,4) = nansum(Alp_mcatch(lid));
    lme_mcatch(L,5) = nansum(Ald_mcatch(lid));
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
    %total biomass
    lme_sbio(L,1) = nansum(Asf_mean(lid));
    lme_sbio(L,2) = nansum(Asp_mean(lid));
    lme_sbio(L,3) = nansum(Asd_mean(lid));
    lme_sbio(L,4) = nansum(Amf_mean(lid));
    lme_sbio(L,5) = nansum(Amp_mean(lid));
    lme_sbio(L,6) = nansum(Amd_mean(lid));
    lme_sbio(L,7) = nansum(Alp_mean(lid));
    lme_sbio(L,8) = nansum(Ald_mean(lid));
    lme_sbio(L,9) = nansum(Ab_mean(lid));
    %total area of LME
    lme_area(L,1) = nansum(AREA_OCN(lid));
end

%% P:D ratio
lme_area_km2 = lme_area * 1e-6;

% POEM LME biomass in MT
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
% MT/km2
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;

clme_mf = NaN*ones(360,200);
clme_mp = clme_mf;
clme_md = clme_mf;
clme_lp = clme_mf;
clme_ld = clme_mf;

plme_AllP = clme_mf;
plme_AllD = clme_mf;

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

    clme_mf(lid) = lme_mcatch(L,1);
    clme_mp(lid) = lme_mcatch(L,2);
    clme_md(lid) = lme_mcatch(L,3);
    clme_lp(lid) = lme_mcatch(L,4);
    clme_ld(lid) = lme_mcatch(L,5);

    lme_sf(lid) = lme_mbio(L,1);
    lme_sp(lid) = lme_mbio(L,2);
    lme_sd(lid) = lme_mbio(L,3);
    lme_mf(lid) = lme_mbio(L,4);
    lme_mp(lid) = lme_mbio(L,5);
    lme_md(lid) = lme_mbio(L,6);
    lme_lp(lid) = lme_mbio(L,7);
    lme_ld(lid) = lme_mbio(L,8);
    lme_b(lid) = lme_mbio(L,9);
    
    plme_AllP(lid) = plme_Pmcatch(L);
    plme_AllD(lid) = plme_Dmcatch(L);
end

clme_AllP = clme_mp+clme_lp;
clme_AllD = clme_md+clme_ld;

lme_AllF = lme_sf+lme_mf;
lme_AllP = lme_sp+lme_mp+lme_lp;
lme_AllD = lme_sd+lme_md+lme_ld;
lme_AllM = lme_mf+lme_mp+lme_md;
lme_AllL = lme_lp+lme_ld;

% Ratios
rPD_biom = lme_AllP ./ (lme_AllP+lme_AllD);
rPF_biom = lme_AllP ./ (lme_AllP+lme_AllF);
rLM_biom = lme_AllL ./ (lme_AllL+lme_AllM);
rPD_catch = clme_AllP ./ (clme_AllP+clme_AllD);
rPD_catch_mtkm2 = plme_AllP ./ (plme_AllP+plme_AllD);


%%
save([dpath 'LME_Hist_1meso_',harv,'_' cfile '.mat'],...
    'lme_mcatch','lme_mbio','lme_sbio','lme_area',...
    'rPD_biom','rPF_biom','rLM_biom','rPD_catch','rPD_catch_mtkm2');



