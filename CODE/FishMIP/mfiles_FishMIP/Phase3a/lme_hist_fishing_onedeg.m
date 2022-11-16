% Calc LME means from output of FEISTY
% 1961-2010 ctrlclim & obsclim with obs fishing effort

clear all
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];

mod = 'obsclim_All_fishobs_v2_';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP/Phase3a/';
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
load([fpath 'Means_Hist_',mod,cfile,'.mat']);

[nid,nt] = size(ld_mean);

%% Map data
%cpath = '/Volumes/MIP/Fish-MIP/Phase3/OneDeg/';
cpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';
load([cpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat']);
load([cpath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat'], 'GRD');
load([cpath 'lme_gfdl-mom6-cobalt2_onedeg.mat'],'tlme');
load([cpath 'cellarea_onedeg.mat'],'cell_area');

[ni,nj]=size(LON);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% Loop over every year to have ts for comp with fishing
close all
lme_mcatch = NaN*ones(66,nt,5);
lme_mbio = NaN*ones(66,nt,9);
lme_sbio = NaN*ones(66,nt,9);

for t=1:nt
    Zsf=NaN*ones(ni,nj);
    Zsp=NaN*ones(ni,nj);
    Zsd=NaN*ones(ni,nj);
    Zmf=NaN*ones(ni,nj);
    Zmp=NaN*ones(ni,nj);
    Zmd=NaN*ones(ni,nj);
    Zlp=NaN*ones(ni,nj);
    Zld=NaN*ones(ni,nj);
    Zb=NaN*ones(ni,nj);
    
    Cmf=NaN*ones(ni,nj);
    Cmp=NaN*ones(ni,nj);
    Cmd=NaN*ones(ni,nj);
    Clp=NaN*ones(ni,nj);
    Cld=NaN*ones(ni,nj);
    
    Zsf(GRD.ID)=sf_mean(:,t);
    Zsp(GRD.ID)=sp_mean(:,t);
    Zsd(GRD.ID)=sd_mean(:,t);
    Zmf(GRD.ID)=mf_mean(:,t);
    Zmp(GRD.ID)=mp_mean(:,t);
    Zmd(GRD.ID)=md_mean(:,t);
    Zlp(GRD.ID)=lp_mean(:,t);
    Zld(GRD.ID)=ld_mean(:,t);
    Zb(GRD.ID)=b_mean(:,t);
    
    Cmf(GRD.ID)=mf_my(:,t);
    Cmp(GRD.ID)=mp_my(:,t);
    Cmd(GRD.ID)=md_my(:,t);
    Clp(GRD.ID)=lp_my(:,t);
    Cld(GRD.ID)=ld_my(:,t);
    
    % g/m2/d --> total g
    Amf_mcatch = Cmf .* cell_area * 365; %mean fish catch per yr
    Amp_mcatch = Cmp .* cell_area * 365;
    Amd_mcatch = Cmd .* cell_area * 365;
    Alp_mcatch = Clp .* cell_area * 365;
    Ald_mcatch = Cld .* cell_area * 365;
    % g/m2 --> total g
    Asf_mean = Zsf .* cell_area;
    Asp_mean = Zsp .* cell_area;
    Asd_mean = Zsd .* cell_area;
    Amf_mean = Zmf .* cell_area;
    Amp_mean = Zmp .* cell_area;
    Amd_mean = Zmd .* cell_area;
    Alp_mean = Zlp .* cell_area;
    Ald_mean = Zld .* cell_area;
    Ab_mean  = Zb .* cell_area;
    
    %% Calc LMEs
    for L=1:66
        lid = find(tlme==L);
        %total catch g
        lme_mcatch(L,t,1) = nansum(Amf_mcatch(lid));
        lme_mcatch(L,t,2) = nansum(Amp_mcatch(lid));
        lme_mcatch(L,t,3) = nansum(Amd_mcatch(lid));
        lme_mcatch(L,t,4) = nansum(Alp_mcatch(lid));
        lme_mcatch(L,t,5) = nansum(Ald_mcatch(lid));
        %mean biomass
        lme_mbio(L,t,1) = nanmean(Asf_mean(lid));
        lme_mbio(L,t,2) = nanmean(Asp_mean(lid));
        lme_mbio(L,t,3) = nanmean(Asd_mean(lid));
        lme_mbio(L,t,4) = nanmean(Amf_mean(lid));
        lme_mbio(L,t,5) = nanmean(Amp_mean(lid));
        lme_mbio(L,t,6) = nanmean(Amd_mean(lid));
        lme_mbio(L,t,7) = nanmean(Alp_mean(lid));
        lme_mbio(L,t,8) = nanmean(Ald_mean(lid));
        lme_mbio(L,t,9) = nanmean(Ab_mean(lid));
        %total biomass
        lme_sbio(L,t,1) = nansum(Asf_mean(lid));
        lme_sbio(L,t,2) = nansum(Asp_mean(lid));
        lme_sbio(L,t,3) = nansum(Asd_mean(lid));
        lme_sbio(L,t,4) = nansum(Amf_mean(lid));
        lme_sbio(L,t,5) = nansum(Amp_mean(lid));
        lme_sbio(L,t,6) = nansum(Amd_mean(lid));
        lme_sbio(L,t,7) = nansum(Alp_mean(lid));
        lme_sbio(L,t,8) = nansum(Ald_mean(lid));
        lme_sbio(L,t,9) = nansum(Ab_mean(lid));
        
    end
end

lme_area = NaN*ones(66,1);
for L=1:66
    lid = find(tlme==L);
    %total area of LME
    lme_area(L,1) = nansum(cell_area(lid));
end

%%
save([fpath 'LME_Hist_',mod,cfile,'.mat'],...
    'lme_mcatch','lme_mbio','lme_sbio','lme_area');
