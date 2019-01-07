% Calc total biomass of POEM
% Climatology
% 150 years
% Saved as mat files

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/v2_kt85_BE75/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);

AREA_OCN = max(area,1);

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

%% Harvested
harv = 'All_fish03';
load([dpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);

% Grid
[ni,nj]=size(lon);

Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zd1=NaN*ones(ni,nj);
Zd2=NaN*ones(ni,nj);
Zm1=NaN*ones(ni,nj);
Zm2=NaN*ones(ni,nj);

Zsf(ID)=sf_mean;
Zsp(ID)=sp_mean;
Zsd(ID)=sd_mean;
Zmf(ID)=mf_mean;
Zmp(ID)=mp_mean;
Zmd(ID)=md_mean;
Zlp(ID)=lp_mean;
Zld(ID)=ld_mean;
Zd1(ID)=all_median1;
Zd2(ID)=all_median2;
Zm1(ID)=all_mean1;
Zm2(ID)=all_mean2;

%g/m2 --> total g
Asf_mean = Zsf .* AREA_OCN;
Asp_mean = Zsp .* AREA_OCN;
Asd_mean = Zsd .* AREA_OCN;
Amf_mean = Zmf .* AREA_OCN;
Amp_mean = Zmp .* AREA_OCN;
Amd_mean = Zmd .* AREA_OCN;
Alp_mean = Zlp .* AREA_OCN;
Ald_mean = Zld .* AREA_OCN;
all_men1 = Zm1 .* AREA_OCN;
all_men2 = Zm2 .* AREA_OCN;
all_med1 = Zd1 .* AREA_OCN;
all_med2 = Zd2 .* AREA_OCN;
 
% Total biomass (compare to J&C15)
tot_bio1 = nansum(all_med1(:)) * 1e-6; %1.5351e+09
tot_bio2 = nansum(all_med2(:)) * 1e-6; %1.5318e+09
tot_bio3 = nansum(all_men1(:)) * 1e-6; %1.5419e+09
tot_bio4 = nansum(all_men2(:)) * 1e-6; %1.5419e+09
tot_bio5 = nansum(Asf_mean(:) + Asp_mean(:) + Asd_mean(:) +...
    Amf_mean(:) + Amp_mean(:) + Amd_mean(:) +...
    Alp_mean(:) + Ald_mean(:)) * 1e-6; %1.5419e+09

tot_bioML = nansum(Amf_mean(:) + Amp_mean(:) + Amd_mean(:) +...
    Alp_mean(:) + Ald_mean(:)) * 1e-6; %1.4983e+09

fish_tot = tot_bio5;
fish_ML = tot_bioML;

%% Pristine
load([dpath 'Means_bio_prod_Climatol_pristine_' cfile '.mat']);

% Grid
Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zd1=NaN*ones(ni,nj);
Zd2=NaN*ones(ni,nj);
Zm1=NaN*ones(ni,nj);
Zm2=NaN*ones(ni,nj);

Zsf(ID)=sf_mean;
Zsp(ID)=sp_mean;
Zsd(ID)=sd_mean;
Zmf(ID)=mf_mean;
Zmp(ID)=mp_mean;
Zmd(ID)=md_mean;
Zlp(ID)=lp_mean;
Zld(ID)=ld_mean;
Zd1(ID)=all_median1;
Zd2(ID)=all_median2;
Zm1(ID)=all_mean1;
Zm2(ID)=all_mean2;

%g/m2 --> total g
Asf_mean = Zsf .* AREA_OCN;
Asp_mean = Zsp .* AREA_OCN;
Asd_mean = Zsd .* AREA_OCN;
Amf_mean = Zmf .* AREA_OCN;
Amp_mean = Zmp .* AREA_OCN;
Amd_mean = Zmd .* AREA_OCN;
Alp_mean = Zlp .* AREA_OCN;
Ald_mean = Zld .* AREA_OCN;
all_men1 = Zm1 .* AREA_OCN;
all_men2 = Zm2 .* AREA_OCN;
all_med1 = Zd1 .* AREA_OCN;
all_med2 = Zd2 .* AREA_OCN;

% Total biomass (compare to J&C15)
unfish_tot = nansum(Asf_mean(:) + Asp_mean(:) + Asd_mean(:) +...
    Amf_mean(:) + Amp_mean(:) + Amd_mean(:) +...
    Alp_mean(:) + Ald_mean(:)) * 1e-6; %1.5419e+09

unfish_ML = nansum(Amf_mean(:) + Amp_mean(:) + Amd_mean(:) +...
    Alp_mean(:) + Ald_mean(:)) * 1e-6; %1.4983e+09

%% Save in table
tab(1,1) = fish_tot;
tab(1,2) = fish_ML;
tab(2,1) = unfish_tot;
tab(2,2) = unfish_ML;
Tab=array2table(tab,'VariableNames',{'All','ML'},...
    'RowNames',{'Fished','Pristine'});
writetable(Tab,[spath 'Tot_biom_pristine_fished_',harv,'_' cfile '.csv'],'Delimiter',',','WriteRowNames',true);
save([dpath 'Tot_biom_pristine_fished_',harv,'_' cfile '.mat'],'tab','Tab');





