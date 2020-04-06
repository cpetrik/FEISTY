% Calc changes from bar plot of trophic amplification
% Ensemble 6, temp3 and orig
% calculate benthic production = benthic biomass * detritus flux * benthic efficiency

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
mzprod_hist = mzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_hist = lzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_hist = npp_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_hist = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
ptemp_hist = ptemp_mean_hist - 273;

mzprod_fore = mzprod_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_fore = lzprod_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_fore = npp_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_fore = det_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
ptemp_fore = ptemp_mean_fore - 273;

zprod_hist = mzprod_hist + lzprod_hist;
zprod_fore = mzprod_fore + lzprod_fore;

% g/m2/d --> total g; Mult flux by 365 - YES
Azprod_hist = 365 * zprod_hist .* AREA_OCN;
Azprod_fore = 365 * zprod_fore .* AREA_OCN;
Adet_hist = 365 * det_hist .* AREA_OCN;
Adet_fore = 365 * det_fore .* AREA_OCN;
Anpp_hist = 365 * npp_hist .* AREA_OCN;
Anpp_fore = 365 * npp_fore .* AREA_OCN;
%nansum(Anpp_mean(:))/9 ~ 50 PgC (1PgC = 1e15gC)


%% Ensemble psets
% epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
%     'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
% load([epath 'simnames_ensem6_mid_kt2_bestAIC_multFup_multPneg.mat'])
nfile = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([nfile 'simnames_ensem6_mid_kt3_bestAIC_multFup_multPneg.mat'])

cF = nan(length(snames)+1 , 1);
cP = nan(length(snames)+1 , 1);
cD = nan(length(snames)+1 , 1);
cS = nan(length(snames)+1 , 1);
cM = nan(length(snames)+1 , 1);
cL = nan(length(snames)+1 , 1);
cB = nan(length(snames)+1 , 1);

hF = nan(length(snames)+1 , 1);
hP = nan(length(snames)+1 , 1);
hD = nan(length(snames)+1 , 1);
hS = nan(length(snames)+1 , 1);
hM = nan(length(snames)+1 , 1);
hL = nan(length(snames)+1 , 1);
hB = nan(length(snames)+1 , 1);

%%
for j = 1:length(snames)
    simname = snames{j};
    fname = pnames{j};
    
    %% Historic Last 50 year means
    load([fname '/Historic_All_fish03_Means_prod_' simname '.mat'],...
        'sf_prod50','sp_prod50','sd_prod50',...
        'mf_prod50','mp_prod50','md_prod50',...
        'lp_prod50','ld_prod50');
    load([fname '/Historic_All_fish03_Means_' simname '.mat'],...
        'b_mean50');
    
    Hsf=sf_prod50;
    Hsp=sp_prod50;
    Hsd=sd_prod50;
    Hmf=mf_prod50;
    Hmp=mp_prod50;
    Hmd=md_prod50;
    Hlp=lp_prod50;
    Hld=ld_prod50;
    % calculate benthic production = benthic biomass * detritus flux * benthic efficiency
    Hb= b_mean50 .* det_hist(ID) .* 0.075;
    
    Hsf(Hsf(:)<0)=0;
    Hsp(Hsp(:)<0)=0;
    Hsd(Hsd(:)<0)=0;
    Hmf(Hmf(:)<0)=0;
    Hmp(Hmp(:)<0)=0;
    Hmd(Hmd(:)<0)=0;
    Hlp(Hlp(:)<0)=0;
    Hld(Hld(:)<0)=0;
    Hb(Hb(:)<0)  =0;
    
    Hsf_area = nansum(Hsf .* AREA_OCN(ID));
    Hsp_area = nansum(Hsp .* AREA_OCN(ID));
    Hsd_area = nansum(Hsd .* AREA_OCN(ID));
    Hmf_area = nansum(Hmf .* AREA_OCN(ID));
    Hmp_area = nansum(Hmp .* AREA_OCN(ID));
    Hmd_area = nansum(Hmd .* AREA_OCN(ID));
    Hlp_area = nansum(Hlp .* AREA_OCN(ID));
    Hld_area = nansum(Hld .* AREA_OCN(ID));
    Hb_area  = nansum(Hb  .* AREA_OCN(ID));
    
    clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50
    clear Hsf Hsp Hsd Hmf Hmp Hmd Hlp Hld Hb
    
    %% Forecast Last 50 year means
    load([fname '/Forecast_All_fish03_Means_prod_' simname '.mat'],...
        'sf_prod50','sp_prod50','sd_prod50',...
        'mf_prod50','mp_prod50','md_prod50',...
        'lp_prod50','ld_prod50');
    load([fname '/Forecast_All_fish03_Means_' simname '.mat'],...
        'b_mean50');
    
    Csf=sf_prod50;
    Csp=sp_prod50;
    Csd=sd_prod50;
    Cmf=mf_prod50;
    Cmp=mp_prod50;
    Cmd=md_prod50;
    Clp=lp_prod50;
    Cld=ld_prod50;
    % calculate benthic production = benthic biomass * detritus flux * benthic efficiency
    Cb= b_mean50 .* det_fore(ID) .* 0.075;
    
    Csf(Csf(:)<0)=0;
    Csp(Csp(:)<0)=0;
    Csd(Csd(:)<0)=0;
    Cmf(Cmf(:)<0)=0;
    Cmp(Cmp(:)<0)=0;
    Cmd(Cmd(:)<0)=0;
    Clp(Clp(:)<0)=0;
    Cld(Cld(:)<0)=0;
    Cb(Cb(:)<0)  =0;
    
    Csf_area = nansum(Csf .* AREA_OCN(ID));
    Csp_area = nansum(Csp .* AREA_OCN(ID));
    Csd_area = nansum(Csd .* AREA_OCN(ID));
    Cmf_area = nansum(Cmf .* AREA_OCN(ID));
    Cmp_area = nansum(Cmp .* AREA_OCN(ID));
    Cmd_area = nansum(Cmd .* AREA_OCN(ID));
    Clp_area = nansum(Clp .* AREA_OCN(ID));
    Cld_area = nansum(Cld .* AREA_OCN(ID));
    Cb_area  = nansum(Cb  .* AREA_OCN(ID));
    
    clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50
    clear Csf Csp Csd Cmf Cmp Cmd Clp Cld Cb
    
    %% Groups
    cF(j) = Csf_area + Cmf_area;
    cP(j) = Csp_area + Cmp_area + Clp_area;
    cD(j) = Csd_area + Cmd_area + Cld_area;
    cS(j) = Csp_area + Csf_area + Csd_area;
    cM(j) = Cmp_area + Cmf_area + Cmd_area;
    cL(j) = Clp_area + Cld_area;
    cB(j) = Cb_area;
    
    hF(j) = Hsf_area + Hmf_area;
    hP(j) = Hsp_area + Hmp_area + Hlp_area ;
    hD(j) = Hsd_area + Hmd_area + Hld_area ;
    hS(j) = Hsp_area + Hsf_area + Hsd_area ;
    hM(j) = Hmp_area + Hmf_area + Hmd_area ;
    hL(j) = Hlp_area + Hld_area;
    hB(j) = Hb_area;
    
end

%% FEISTY Output orig pset
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];

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
Hb = b_mean50 .* det_hist(ID) .* 0.075;

Hsf(Hsf(:)<0)=0;
Hsp(Hsp(:)<0)=0;
Hsd(Hsd(:)<0)=0;
Hmf(Hmf(:)<0)=0;
Hmp(Hmp(:)<0)=0;
Hmd(Hmd(:)<0)=0;
Hlp(Hlp(:)<0)=0;
Hld(Hld(:)<0)=0;
Hb(Hb(:)<0)  =0;

Hsf_area = nansum(Hsf .* AREA_OCN(ID));
Hsp_area = nansum(Hsp .* AREA_OCN(ID));
Hsd_area = nansum(Hsd .* AREA_OCN(ID));
Hmf_area = nansum(Hmf .* AREA_OCN(ID));
Hmp_area = nansum(Hmp .* AREA_OCN(ID));
Hmd_area = nansum(Hmd .* AREA_OCN(ID));
Hlp_area = nansum(Hlp .* AREA_OCN(ID));
Hld_area = nansum(Hld .* AREA_OCN(ID));
Hb_area  = nansum(Hb  .* AREA_OCN(ID));

clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50
clear Hsf Hsp Hsd Hmf Hmp Hmd Hlp Hld Hb
    
% Forecast
load([fpath 'Means_fore_',harv,'_' cfile '.mat'],...
    'sf_prod50','sp_prod50','sd_prod50',...
    'mf_prod50','mp_prod50','md_prod50',...
    'lp_prod50','ld_prod50','b_mean50');

Csf=sf_prod50;
Csp=sp_prod50;
Csd=sd_prod50;
Cmf=mf_prod50;
Cmp=mp_prod50;
Cmd=md_prod50;
Clp=lp_prod50;
Cld=ld_prod50;
Cb = b_mean50 .* det_fore(ID) .* 0.075;

Csf(Csf(:)<0)=0;
Csp(Csp(:)<0)=0;
Csd(Csd(:)<0)=0;
Cmf(Cmf(:)<0)=0;
Cmp(Cmp(:)<0)=0;
Cmd(Cmd(:)<0)=0;
Clp(Clp(:)<0)=0;
Cld(Cld(:)<0)=0;
Cb(Cb(:)<0)  =0;

Csf_area = nansum(Csf .* AREA_OCN(ID));
Csp_area = nansum(Csp .* AREA_OCN(ID));
Csd_area = nansum(Csd .* AREA_OCN(ID));
Cmf_area = nansum(Cmf .* AREA_OCN(ID));
Cmp_area = nansum(Cmp .* AREA_OCN(ID));
Cmd_area = nansum(Cmd .* AREA_OCN(ID));
Clp_area = nansum(Clp .* AREA_OCN(ID));
Cld_area = nansum(Cld .* AREA_OCN(ID));
Cb_area  = nansum(Cb  .* AREA_OCN(ID));

clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50
clear Csf Csp Csd Cmf Cmp Cmd Clp Cld Cb
    
% Groups
np=length(snames)+1;
cF(np) = Csf_area + Cmf_area;
cP(np) = Csp_area + Cmp_area + Clp_area;
cD(np) = Csd_area + Cmd_area + Cld_area;
cS(np) = Csp_area + Csf_area + Csd_area;
cM(np) = Cmp_area + Cmf_area + Cmd_area;
cL(np) = Clp_area + Cld_area;
cB(np) = Cb_area;

hF(np) = Hsf_area + Hmf_area;
hP(np) = Hsp_area + Hmp_area + Hlp_area ;
hD(np) = Hsd_area + Hmd_area + Hld_area ;
hS(np) = Hsp_area + Hsf_area + Hsd_area ;
hM(np) = Hmp_area + Hmf_area + Hmd_area ;
hL(np) = Hlp_area + Hld_area;
hB(np) = Hb_area;

%% Percent differences
hAll = hF+hP+hD;
hPel = hF+hP;
cAll = cF+cP+cD;
cPel = cF+cP;

pbar(1) = (nansum(Anpp_fore(:))-nansum(Anpp_hist(:))) ./ nansum(Anpp_hist(:));
pbar(2) = (nansum(Azprod_fore(:))-nansum(Azprod_hist(:))) ./ nansum(Azprod_hist(:)); %-0.079
pbar(3) = (nansum(Adet_fore(:))-nansum(Adet_hist(:))) ./ nansum(Adet_hist(:));
pnames = {'NPP','Z','Det'};

fbar(:,1) = ((cM(:))-(hM(:))) ./ (hM(:));
fbar(:,2) = ((cL(:))-(hL(:))) ./ (hL(:));
fbar(:,3) = ((cF(:))-(hF(:))) ./ (hF(:));
fbar(:,4) = ((cP(:))-(hP(:))) ./ (hP(:));
fbar(:,5) = ((cD(:))-(hD(:))) ./ (hD(:));
fbar(:,6) = ((cPel(:))-(hPel(:))) ./ (hPel(:));
fbar(:,7) = ((cAll(:))-(hAll(:))) ./ (hAll(:));
fbar(:,8) = ((cB(:))-(hB(:))) ./ (hB(:));
names = {'M','L','F','P','D','Pel','All','B'};

save([nfile 'Prod_diff_50yr_ensem6_mid_kt3_bestAIC_multFup_multPneg.mat'],...
    'hF','hP','hD','hS','hM','hL','hAll','hPel','hB',...
    'cF','cP','cD','cS','cM','cL','cAll','cPel','cB',...
    'pbar','fbar','pnames','names');


