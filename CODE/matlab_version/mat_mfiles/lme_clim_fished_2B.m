% Calc LME biomass of POEM
% Climatology
% 150 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/FEISTY/NC/Matlab_new_size/';

Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
cdir='/Volumes/FEISTY/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cdir 'temp_100_1deg_ESM26_5yr_clim_191_195.mat'])
load([cdir 'btm_temp_1deg_ESM26_5yr_clim_191_195.mat'])

ptemp_mean_clim=squeeze(nanmean(temp_100,1));
btemp_mean_clim=squeeze(nanmean(btm_temp,1));

%%
AREA_OCN = max(area,1);

%cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_2B_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

%load([dpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);
load([dpath 'Means_Climatol_' harv '_' cfile '.mat']);

%% Plots in space
[ni,nj]=size(lon);

Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zsb=NaN*ones(ni,nj);
Zmb=NaN*ones(ni,nj);

Tsf=NaN*ones(ni,nj);
Tsp=NaN*ones(ni,nj);
Tsd=NaN*ones(ni,nj);
Tmf=NaN*ones(ni,nj);
Tmp=NaN*ones(ni,nj);
Tmd=NaN*ones(ni,nj);
Tlp=NaN*ones(ni,nj);
Tld=NaN*ones(ni,nj);

Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);

Zsf(ID)=sf_mean;
Zsp(ID)=sp_mean;
Zsd(ID)=sd_mean;
Zmf(ID)=mf_mean;
Zmp(ID)=mp_mean;
Zmd(ID)=md_mean;
Zlp(ID)=lp_mean;
Zld(ID)=ld_mean;
Zsb(ID)=sb_mean;
Zmb(ID)=mb_mean;

Tsf(ID)=sf_tot;
Tsp(ID)=sp_tot;
Tsd(ID)=sd_tot;
Tmf(ID)=mf_tot;
Tmp(ID)=mp_tot;
Tmd(ID)=md_tot;
Tlp(ID)=lp_tot;
Tld(ID)=ld_tot;

Cmf(ID)=mf_my;
Cmp(ID)=mp_my;
Cmd(ID)=md_my;
Clp(ID)=lp_my;
Cld(ID)=ld_my;

% g/m2/d --> total g
Amf_mcatch = Cmf .* AREA_OCN * 365; %mean fish catch per yr
Amp_mcatch = Cmp .* AREA_OCN * 365;
Amd_mcatch = Cmd .* AREA_OCN * 365;
Alp_mcatch = Clp .* AREA_OCN * 365;
Ald_mcatch = Cld .* AREA_OCN * 365;
% g/m2 --> total g
Asf_mean = Zsf .* AREA_OCN;
Asp_mean = Zsp .* AREA_OCN;
Asd_mean = Zsd .* AREA_OCN;
Amf_mean = Zmf .* AREA_OCN;
Amp_mean = Zmp .* AREA_OCN;
Amd_mean = Zmd .* AREA_OCN;
Alp_mean = Zlp .* AREA_OCN;
Ald_mean = Zld .* AREA_OCN;
Asb_mean = Zsb .* AREA_OCN;
Amb_mean = Zmb .* AREA_OCN;

Asf_tot = Tsf .* AREA_OCN;
Asp_tot = Tsp .* AREA_OCN;
Asd_tot = Tsd .* AREA_OCN;
Amf_tot = Tmf .* AREA_OCN;
Amp_tot = Tmp .* AREA_OCN;
Amd_tot = Tmd .* AREA_OCN;
Alp_tot = Tlp .* AREA_OCN;
Ald_tot = Tld .* AREA_OCN;

%% Calc LMEs
tlme = lme_mask_onedeg;

lme_mcatch = NaN*ones(66,5);
lme_mbio = NaN*ones(66,10);
lme_sbio = NaN*ones(66,10);
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
    lme_mbio(L,9) = nanmean(Asb_mean(lid));
    lme_mbio(L,10) = nanmean(Amb_mean(lid));
    
    %total biomass
    lme_sbio(L,1) = nansum(Asf_mean(lid));
    lme_sbio(L,2) = nansum(Asp_mean(lid));
    lme_sbio(L,3) = nansum(Asd_mean(lid));
    lme_sbio(L,4) = nansum(Amf_mean(lid));
    lme_sbio(L,5) = nansum(Amp_mean(lid));
    lme_sbio(L,6) = nansum(Amd_mean(lid));
    lme_sbio(L,7) = nansum(Alp_mean(lid));
    lme_sbio(L,8) = nansum(Ald_mean(lid));
    lme_sbio(L,9) = nansum(Asb_mean(lid));
    lme_sbio(L,10) = nansum(Amb_mean(lid));
    %total area of LME
    lme_area(L,1) = nansum(AREA_OCN(lid));
end

%catch in lmes
Tlme_mcatch = sum(lme_mcatch(:));
tot_catch = (Amf_mcatch+Amp_mcatch+Amd_mcatch+Alp_mcatch+Ald_mcatch);
tot_catch2 = nansum(tot_catch(:));
frac_mcatch_lme = Tlme_mcatch/tot_catch2

%%
save([dpath 'LME_clim_fished_',harv,'_' cfile '.mat'],...
    'lme_mcatch','lme_mbio','lme_sbio','lme_area');

