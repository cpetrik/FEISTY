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

exper = 'CORE_Hindcast1988_no_move_All_fish03_yield';
load([dpath 'Annual_Means_' exper '_' cfile '.mat']);

%% Total annual catch - reshape to grid
% e.g.  MF.catch = MF.yield .*mos .*area_mat;
%       units_catch = 'g_mo';
%       mf_tac = sum(MF.catch(:,st(n):en(n)),2,'omitnan');
%       units_tac = 'g_yr';

[ni,nj]  = size(GRD2.area);
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
vlme = tlme(GRD1.ID);

lme_tcatch_MF = nan*ones(66,nt);
lme_tcatch_MP = lme_tcatch_MF;
lme_tcatch_MD = lme_tcatch_MF;
lme_tcatch_LP = lme_tcatch_MF;
lme_tcatch_LD = lme_tcatch_MF;

for L=1:66
    lid = find(vlme==L);
    %total catch g
    lme_tcatch_MF(L,:) = sum(mf_tac(lid,:),'omitnan');
    lme_tcatch_MP(L,:) = sum(mp_tac(lid,:),'omitnan');
    lme_tcatch_MD(L,:) = sum(md_tac(lid,:),'omitnan');
    lme_tcatch_LP(L,:) = sum(lp_tac(lid,:),'omitnan');
    lme_tcatch_LD(L,:) = sum(ld_tac(lid,:),'omitnan');
end

lme_tcatch_all = lme_tcatch_MF + lme_tcatch_MP + lme_tcatch_MD +...
    lme_tcatch_LP + lme_tcatch_LD;
lme_tcatch_MT = sum(lme_tcatch_all,'omitnan') * 1e-6; %g to MT?

%%
save([dpath 'LME_Annual_total_catch_' exper '_' cfile '.mat'],...
    'lme_tcatch_MF','lme_tcatch_MP','lme_tcatch_MD',...
    'lme_tcatch_LP','lme_tcatch_LD','lme_tcatch_all','lme_area');





