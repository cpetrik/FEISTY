% save P:D ratio on grid
% Climatology
% 150 years
% Saved as mat files

% Calc LME biomass of POEM
% Climatology
% 150 years
% Saved as mat files

clear all
close all

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

ptemp_mean_clim=squeeze(nanmean(temp_100,1));
btemp_mean_clim=squeeze(nanmean(btm_temp,1));

%%
AREA_OCN = max(area,1);

%cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);

%% Plots in space
[ni,nj]=size(lon);

Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);

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

F_catch_gm2d = Cmf;
P_catch_gm2d = Cmp+Clp;
D_catch_gm2d = Cmd+Cld;

F_catch_g = Amf_mcatch;
P_catch_g = Amp_mcatch+Alp_mcatch;
D_catch_g = Amd_mcatch+Ald_mcatch;


%%
save([dpath 'PD_Clim_fished_',harv,'_' cfile '.mat'],...
    'F_catch_gm2d','P_catch_gm2d','D_catch_gm2d',...
    'F_catch_g','P_catch_g','D_catch_g',...
    'lat','lon','AREA_OCN');

%% csv
vec(:,1)=lat(:);
vec(:,2)=lon(:);
vec(:,3)=AREA_OCN(:);
vec(:,4)=F_catch_gm2d(:);
vec(:,5)=P_catch_gm2d(:);
vec(:,6)=D_catch_gm2d(:);
vec(:,7)=F_catch_g(:);
vec(:,8)=P_catch_g(:);
vec(:,9)=D_catch_g(:);

oc=find(vec(:,3)~=1);
vec = vec(oc,:);

VecT = array2table(vec,'VariableNames',{'lat','lon','area',...
    'F_catch_gm2d','P_catch_gm2d','D_catch_gm2d',...
    'F_catch_g','P_catch_g','D_catch_g'});
writetable(VecT,[dpath 'PD_Clim_fished_',harv,'_' cfile '.csv'],'Delimiter',',','WriteRowNames',true)

