% Calc LME biomass of FEISTY
% CORE-forced Hindcast 
% Saved as mat files

clear 
close all

%%
vpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CORE-forced/';
cpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/grid_cobalt/';

%1-D
load([vpath 'Data_grid_ocean_cobalt_ESM2Mcore.mat'],'GRD');
GRD1 = GRD;
clear GRD

%2-D
load([vpath 'Data_hindcast_grid_cp2D.mat'],'GRD')
GRD2 = GRD;
clear GRD

%%
load([vpath 'lme_mask_esm2m.mat']);


%% FEISTY
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/';
dp = '/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_CC80_RE00100';

dpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/CORE/'];
%dpath=['/project/Feisty/NC/Matlab_new_size/' cfile '/CORE/'];

ppath = [pp cfile '/CORE/'];

load([dpath 'Annual_Means_' exper '_' cfile '.mat']);

%% Total annual catch - reshape to grid
% e.g.  MF.catch = MF.yield .*mos .*area_mat;
%       units_catch = 'g_mo';
%       mf_tac = sum(MF.catch(:,st(n):en(n)),2,'omitnan');
%       units_tac = 'g_yr';

[ni,nj]=size(geolon_t);
[nid,nt]=size(mf_tac);

%% Calc LMEs
tlme = lme_mask_esm2m';

lme_area = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    %total area of LME
    lme_area(L,1) = sum(GRD2.area(lid),'omitnan');
end

%% g/yr, area already accounted for
% can just sum over lme
vlme = tlme(ID);

lme_tcatch_MF = nan*ones(66,nyr);
lme_tcatch_MP = lme_tcatch_MF;
lme_tcatch_MD = lme_tcatch_MF;
lme_tcatch_LP = lme_tcatch_MF;
lme_tcatch_LD = lme_tcatch_MF;

for L=1:66
    lid = find(vlme==L);
    %total catch g
    lme_tcatch_MF(L,:) = sum(Amf_mean(lid,:),'omitnan');
    lme_tcatch_MP(L,:) = sum(Amp_mean(lid,:),'omitnan');
    lme_tcatch_MD(L,:) = sum(Amd_mean(lid,:),'omitnan');
    lme_tcatch_LP(L,:) = sum(Alp_mean(lid,:),'omitnan');
    lme_tcatch_LD(L,:) = sum(Ald_mean(lid,:),'omitnan');
end

%hlme_tcatch = sum(lme_tcatch,2,'omitnan') * 1e-6; %g to MT?

%%
save([dpath 'LME_Annual_total_catch_' exper '_' cfile '.mat'],...
    'lme_tcatch_MF','lme_tcatch_MP','lme_tcatch_MD',...
    'lme_tcatch_LP','lme_tcatch_LD','lme_area');





