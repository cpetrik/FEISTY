% Calc changes from bar plot of trophic amplification
% Ensemble and orig

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

load([cpath 'hindcast_gridspec.mat'],'AREA_OCN');
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

%% NPP and zoop
load([bpath 'cobalt_det_temp_zoop_npp_means.mat']);

% molN/m2/s --> g/m2/d
mzloss_hist = mzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_hist = lzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_hist = npp_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_hist = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
ptemp_hist = ptemp_mean_hist - 273;

mzloss_fore = mzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_fore = lzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_fore = npp_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_fore = det_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
ptemp_fore = ptemp_mean_fore - 273;

zloss_hist = mzloss_hist + lzloss_hist;
zloss_fore = mzloss_fore + lzloss_fore;

% g/m2/d --> total g; Mult flux by 365 - YES
Azloss_hist = 365 * zloss_hist .* AREA_OCN;
Azloss_fore = 365 * zloss_fore .* AREA_OCN;
Adet_hist = 365 * det_hist .* AREA_OCN;
Adet_fore = 365 * det_fore .* AREA_OCN;
Anpp_hist = 365 * npp_hist .* AREA_OCN;
Anpp_fore = 365 * npp_fore .* AREA_OCN;
%nansum(Anpp_mean(:))/9 ~ 50 PgC (1PgC = 1e15gC)


%% Ensemble psets
epath = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([epath 'simnames_ensem5_mid5_bestAIC_multFup_multPneg.mat'])

cF = nan(length(ID),length(snames)+1);
cP = nan(length(ID),length(snames)+1);
cD = nan(length(ID),length(snames)+1);
cS = nan(length(ID),length(snames)+1);
cM = nan(length(ID),length(snames)+1);
cL = nan(length(ID),length(snames)+1);

hF = nan(length(ID),length(snames)+1);
hP = nan(length(ID),length(snames)+1);
hD = nan(length(ID),length(snames)+1);
hS = nan(length(ID),length(snames)+1);
hM = nan(length(ID),length(snames)+1);
hL = nan(length(ID),length(snames)+1);

%%
for j = 1:length(snames)
    simname = snames{j};
    fname = pnames{j};
    
    %% Historic Last 50 year means
    load([fname '/Historic_All_fish03_Means_prod_' simname '.mat'],...
        'sf_prod50','sp_prod50','sd_prod50',...
        'mf_prod50','mp_prod50','md_prod50',...
        'lp_prod50','ld_prod50');
    
    Hsf=sf_prod50;
    Hsp=sp_prod50;
    Hsd=sd_prod50;
    Hmf=mf_prod50;
    Hmp=mp_prod50;
    Hmd=md_prod50;
    Hlp=lp_prod50;
    Hld=ld_prod50;
    
    Hsf(Hsf(:)<0)=0;
    Hsp(Hsp(:)<0)=0;
    Hsd(Hsd(:)<0)=0;
    Hmf(Hmf(:)<0)=0;
    Hmp(Hmp(:)<0)=0;
    Hmd(Hmd(:)<0)=0;
    Hlp(Hlp(:)<0)=0;
    Hld(Hld(:)<0)=0;
    
    Hsf_area = Hsf .* AREA_OCN(ID);
    Hsp_area = Hsp .* AREA_OCN(ID);
    Hsd_area = Hsd .* AREA_OCN(ID);
    Hmf_area = Hmf .* AREA_OCN(ID);
    Hmp_area = Hmp .* AREA_OCN(ID);
    Hmd_area = Hmd .* AREA_OCN(ID);
    Hlp_area = Hlp .* AREA_OCN(ID);
    Hld_area = Hld .* AREA_OCN(ID);
    
    clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50
    clear Hsf Hsp Hsd Hmf Hmp Hmd Hlp Hld
    
    %% Forecast Last 50 year means
    load([fname '/Forecast_All_fish03_Means_prod_' simname '.mat'],...
        'sf_prod50','sp_prod50','sd_prod50',...
        'mf_prod50','mp_prod50','md_prod50',...
        'lp_prod50','ld_prod50');
    
    Csf(ID)=sf_prod50;
    Csp(ID)=sp_prod50;
    Csd(ID)=sd_prod50;
    Cmf(ID)=mf_prod50;
    Cmp(ID)=mp_prod50;
    Cmd(ID)=md_prod50;
    Clp(ID)=lp_prod50;
    Cld(ID)=ld_prod50;
    
    Csf(Csf(:)<0)=0;
    Csp(Csp(:)<0)=0;
    Csd(Csd(:)<0)=0;
    Cmf(Cmf(:)<0)=0;
    Cmp(Cmp(:)<0)=0;
    Cmd(Cmd(:)<0)=0;
    Clp(Clp(:)<0)=0;
    Cld(Cld(:)<0)=0;
    
    Csf_area = Csf .* AREA_OCN(ID);
    Csp_area = Csp .* AREA_OCN(ID);
    Csd_area = Csd .* AREA_OCN(ID);
    Cmf_area = Cmf .* AREA_OCN(ID);
    Cmp_area = Cmp .* AREA_OCN(ID);
    Cmd_area = Cmd .* AREA_OCN(ID);
    Clp_area = Clp .* AREA_OCN(ID);
    Cld_area = Cld .* AREA_OCN(ID);
    
    clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50
    clear Csf Csp Csd Cmf Cmp Cmd Clp Cld
    
    %% Groups
    cF(:,j) = Csf_area + Cmf_area;
    cP(:,j) = Csp_area + Cmp_area + Clp_area;
    cD(:,j) = Csd_area + Cmd_area + Cld_area;
    cS(:,j) = Csp_area + Csf_area + Csd_area;
    cM(:,j) = Cmp_area + Cmf_area + Cmd_area;
    cL(:,j) = Clp_area + Cld_area;
    
    hF(:,j) = Hsf_area + Hmf_area;
    hP(:,j) = Hsp_area + Hmp_area + Hlp_area ;
    hD(:,j) = Hsd_area + Hmd_area + Hld_area ;
    hS(:,j) = Hsp_area + Hsf_area + Hsd_area ;
    hM(:,j) = Hmp_area + Hmf_area + Hmd_area ;
    hL(:,j) = Hlp_area + Hld_area;
    
end

%% FEISTY Output orig pset
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];

% Hindcast
load([fpath 'Means_Historic_' harv '_prod_' cfile '.mat'],...
    'sf_prod50','sp_prod50','sd_prod50',...
    'mf_prod50','mp_prod50','md_prod50',...
    'lp_prod50','ld_prod50','b_mean50');

Hsf=sf_prod50;
Hsp=sp_prod50;
Hsd=sd_prod50;
Hmf=mf_prod50;
Hmp=mp_prod50;
Hmd=md_prod50;
Hlp=lp_prod50;
Hld=ld_prod50;

Hsf(Hsf(:)<0)=0;
Hsp(Hsp(:)<0)=0;
Hsd(Hsd(:)<0)=0;
Hmf(Hmf(:)<0)=0;
Hmp(Hmp(:)<0)=0;
Hmd(Hmd(:)<0)=0;
Hlp(Hlp(:)<0)=0;
Hld(Hld(:)<0)=0;

Hsf_area = Hsf .* AREA_OCN(ID);
Hsp_area = Hsp .* AREA_OCN(ID);
Hsd_area = Hsd .* AREA_OCN(ID);
Hmf_area = Hmf .* AREA_OCN(ID);
Hmp_area = Hmp .* AREA_OCN(ID);
Hmd_area = Hmd .* AREA_OCN(ID);
Hlp_area = Hlp .* AREA_OCN(ID);
Hld_area = Hld .* AREA_OCN(ID);

clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50
clear Hsf Hsp Hsd Hmf Hmp Hmd Hlp Hld
    
% Forecast
load([fpath 'Means_fore_',harv,'_' cfile '.mat'],...
    'sf_prod50','sp_prod50','sd_prod50',...
    'mf_prod50','mp_prod50','md_prod50',...
    'lp_prod50','ld_prod50','b_mean50');

Csf(ID)=sf_prod50;
Csp(ID)=sp_prod50;
Csd(ID)=sd_prod50;
Cmf(ID)=mf_prod50;
Cmp(ID)=mp_prod50;
Cmd(ID)=md_prod50;
Clp(ID)=lp_prod50;
Cld(ID)=ld_prod50;

Csf(Csf(:)<0)=0;
Csp(Csp(:)<0)=0;
Csd(Csd(:)<0)=0;
Cmf(Cmf(:)<0)=0;
Cmp(Cmp(:)<0)=0;
Cmd(Cmd(:)<0)=0;
Clp(Clp(:)<0)=0;
Cld(Cld(:)<0)=0;

Csf_area = Csf .* AREA_OCN;
Csp_area = Csp .* AREA_OCN;
Csd_area = Csd .* AREA_OCN;
Cmf_area = Cmf .* AREA_OCN;
Cmp_area = Cmp .* AREA_OCN;
Cmd_area = Cmd .* AREA_OCN;
Clp_area = Clp .* AREA_OCN;
Cld_area = Cld .* AREA_OCN;

clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50
clear Csf Csp Csd Cmf Cmp Cmd Clp Cld
    
% Groups
cF(:,15) = Csf_area + Cmf_area;
cP(:,15) = Csp_area + Cmp_area + Clp_area;
cD(:,15) = Csd_area + Cmd_area + Cld_area;
cS(:,15) = Csp_area + Csf_area + Csd_area;
cM(:,15) = Cmp_area + Cmf_area + Cmd_area;
cL(:,15) = Clp_area + Cld_area;

hF(:,15) = Hsf_area + Hmf_area;
hP(:,15) = Hsp_area + Hmp_area + Hlp_area ;
hD(:,15) = Hsd_area + Hmd_area + Hld_area ;
hS(:,15) = Hsp_area + Hsf_area + Hsd_area ;
hM(:,15) = Hmp_area + Hmf_area + Hmd_area ;
hL(:,15) = Hlp_area + Hld_area;

%% Percent differences

pbar(1) = (nansum(Anpp_fore(:))-nansum(Anpp_hist(:))) ./ nansum(Anpp_hist(:));
pbar(2) = -0.079; %(nansum(Azloss_fore(:))-nansum(Azloss_hist(:))) ./ nansum(Azloss_hist(:));
pbar(3) = (nansum(Adet_fore(:))-nansum(Adet_hist(:))) ./ nansum(Adet_hist(:));
pbar(4) = (nansum(cM(:))-nansum(hM(:))) ./ nansum(hM(:));
pbar(5) = (nansum(cL(:))-nansum(hL(:))) ./ nansum(hL(:));
pbar(6) = (nansum(cF(:))-nansum(hF(:))) ./ nansum(hF(:));
pbar(7) = (nansum(cP(:))-nansum(hP(:))) ./ nansum(hP(:));
pbar(8) = (nansum(cD(:))-nansum(hD(:))) ./ nansum(hD(:));
pbar(9) = (nansum(cPel(:))-nansum(hPel(:))) ./ nansum(hPel(:));
pbar(10) = (nansum(cAll(:))-nansum(hAll(:))) ./ nansum(hAll(:));
pnames = {'NPP','Z','Det','M','L','F','P','D','Pel','All'};

delB = (nansum(Cb_area(:))-nansum(Hb_area(:))) ./ nansum(Hb_area(:));
%change in B is -10.16%, but change in Det is -12.3%

