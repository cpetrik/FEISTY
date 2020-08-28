function [misf] = hist_1meso_ratios_zoo_misfit_loop(sF,sP,sD,sB,sM,sL,sZ)

%1 meso FEISTY w/o overcon vs. w/overcon
%ratios of P:D, P:F, L:M
%frac of HPloss consumed

%% 1 meso structurally FEISTY
cfile1 = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
fpath1 = ['/Volumes/FEISTY/NC/Matlab_new_size/' cfile1 '/'];
harv = 'All_fish03';

load([fpath1 'Means_Historic_1meso_',harv,'_' cfile1 '.mat'],...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50','b_mean50','mz_mfrac5');

CF = sf_mean50 + mf_mean50;
CP = sp_mean50 + mp_mean50 + lp_mean50;
CD = sd_mean50 + md_mean50 + ld_mean50;
CB = b_mean50;
CM = mf_mean50 + mp_mean50 + md_mean50;
CL = lp_mean50 + ld_mean50;
CZ = mz_mfrac5;

CAll = CF+CP+CD;
cFracPD = CP ./ (CP+CD);
cFracPF = CP ./ (CP+CF);
cFracLM = CL ./ (CL+CM);

%% 1 meso without overcon
ZF=sF;
ZP=sP;
ZD=sD;
ZM=sM;
ZL=sL;
ZB=sB;
ZZ=sZ;

ZAll = ZF+ZP+ZD;
zFracPD = ZP ./ (ZP+ZD);
zFracPF = ZP ./ (ZP+ZF);
zFracLM = ZL ./ (ZL+ZM);

%% Stats -------------------------
misf = nan*ones(length(sF),9);

misf(:,1) = (ZF-CF);
misf(:,2) = (ZP-CP);
misf(:,3) = (ZD-CD);
misf(:,5) = (ZB-CB);
misf(:,4) = (ZAll-CAll);
misf(:,6) = zFracPD - cFracPD;
misf(:,7) = zFracPF - cFracPF;
misf(:,8) = zFracLM - cFracLM;
misf(:,9) = (ZZ-CZ);


