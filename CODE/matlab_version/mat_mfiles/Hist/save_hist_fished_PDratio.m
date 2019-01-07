% Calc LME biomass of POEM
% Hindcast time period 1990-1995 (for comparing with Climatology)
% 5 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t','AREA_OCN');
grid = csvread([cpath 'grid_csv.csv']);
load([cpath 'lme_mask_esm2m.mat']);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

%% FEISTY
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'Means_Historic_',harv,'_' cfile '.mat']);

%% Plots in space
[ni,nj]=size(geolon_t);

Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);

Zmf(grid(:,1))=mf_my5;
Zmp(grid(:,1))=mp_my5;
Zmd(grid(:,1))=md_my5;
Zlp(grid(:,1))=lp_my5;
Zld(grid(:,1))=ld_my5;

% g/m2 --> total g
Amf_mcatch = Zmf .* AREA_OCN * 365; %mean fish catch per yr
Amp_mcatch = Zmp .* AREA_OCN * 365;
Amd_mcatch = Zmd .* AREA_OCN * 365;
Alp_mcatch = Zlp .* AREA_OCN * 365;
Ald_mcatch= Zld .* AREA_OCN * 365;

F_catch_gm2d = Zmf;
P_catch_gm2d = Zmp+Zlp;
D_catch_gm2d = Zmd+Zld;

F_catch_g = Amf_mcatch;
P_catch_g = Amp_mcatch+Alp_mcatch;
D_catch_g = Amd_mcatch+Ald_mcatch;


%%
save([dpath 'PD_Clim_fished_',harv,'_' cfile '.mat'],...
    'F_catch_gm2d','P_catch_gm2d','D_catch_gm2d',...
    'F_catch_g','P_catch_g','D_catch_g',...
    'geolat_t','geolon_t','AREA_OCN');

%% csv
vec(:,1)=geolat_t(:);
vec(:,2)=geolon_t(:);
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

