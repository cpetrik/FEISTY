% Check that new version with param structure gives same results

clear all
close all

load('/Volumes/FEISTY/FEISTY_clim/ESM26_1deg_5yr_clim_191_195_grid.mat','GRD');
NX = length(GRD.Z);
ID = 1:NX;
param.NX = NX;
param.ID = ID;

load('/Volumes/FEISTY/FEISTY_clim/LHS_param6_mid6_samek.mat','fx');

% PARAMETER SENSITIVITY TEST
j = 1;

%! Change individual parameters
pset = fx(j,:);
%pset=round(pset,3);
param = set_params6_samek(pset,param);

%! Make core parameters/constants (global)
param = const_params6_samek(param);

%! Create a directory for output
fname = sub_fname_ensemble6_samek(param);

%% Last year means
load([fname,'.mat'])

[id,nt] = size(Clim_Bent.bio);
time=1:nt;
lyr=time((end-12+1):end);
sp_mean1=mean(Clim_Sml_p.bio(:,lyr),2);
sf_mean1=mean(Clim_Sml_f.bio(:,lyr),2);
sd_mean1=mean(Clim_Sml_d.bio(:,lyr),2);
mp_mean1=mean(Clim_Med_p.bio(:,lyr),2);
mf_mean1=mean(Clim_Med_f.bio(:,lyr),2);
md_mean1=mean(Clim_Med_d.bio(:,lyr),2);
lp_mean1=mean(Clim_Lrg_p.bio(:,lyr),2);
ld_mean1=mean(Clim_Lrg_d.bio(:,lyr),2);
b_mean1=mean(Clim_Bent.bio(:,lyr),2);

%%
fname2 = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/Climatol_All_fish030_Lam600_enc50-b150-k063_cmax-k063_met300-b150-k063';
load([fname2,'.mat'])

[id,nt] = size(Clim_Bent.bio);
time=1:nt;
lyr=time((end-12+1):end);
sp_mean2=mean(Clim_Sml_p.bio(:,lyr),2);
sf_mean2=mean(Clim_Sml_f.bio(:,lyr),2);
sd_mean2=mean(Clim_Sml_d.bio(:,lyr),2);
mp_mean2=mean(Clim_Med_p.bio(:,lyr),2);
mf_mean2=mean(Clim_Med_f.bio(:,lyr),2);
md_mean2=mean(Clim_Med_d.bio(:,lyr),2);
lp_mean2=mean(Clim_Lrg_p.bio(:,lyr),2);
ld_mean2=mean(Clim_Lrg_d.bio(:,lyr),2);
b_mean2=mean(Clim_Bent.bio(:,lyr),2);

%%
sfs = sum(sf_mean1==sf_mean2);
sps = sum(sp_mean1==sp_mean2);
sds = sum(sd_mean1==sd_mean2);
mfs = sum(mf_mean1==mf_mean2);
mps = sum(mp_mean1==mp_mean2);
mds = sum(md_mean1==md_mean2);
lps = sum(lp_mean1==lp_mean2);
lds = sum(ld_mean1==ld_mean2);
bs = sum(b_mean1==b_mean2);

