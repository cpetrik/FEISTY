% Get fishing yield results from Terra runs

clear all
close all

load('ESM26_1deg_5yr_clim_191_195_grid.mat','GRD');
load('ESM26_1deg_5yr_clim_191_195_gridspec.mat','ID');

area_ocn=load('esm26_area_1deg.mat');
AREA_OCN = max(area_ocn.area,1);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);
area = AREA_OCN(ID);
area_km2 = area * 1e-6;

%% List all outputs currently available
fpath = '/scratch/user/cpetrik/FEISTY_output/Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/';

Files = dir([fpath 'Climatol_fish_F*.mat']);
Fnames = {Files.name};
Fpaths = {Files.folder};

%% Last year means

MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

mf_tc = NaN*ones(length(ID),length(Fnames));
mp_tc = mf_tc;
md_tc = mf_tc;
lp_tc = mf_tc;
ld_tc = mf_tc;
mf_gtc = NaN*ones(length(Fnames),1);
mp_gtc = mf_gtc;
md_gtc = mf_gtc;
lp_gtc = mf_gtc;
ld_gtc = mf_gtc;

for n = 1:5 %1:length(Fnames)
    load([Fpaths{n} '/' Fnames{n}]);
    
    [id,nt] = size(Clim_Lrg_d.yield);
    time=1:nt;
    lyr=time((end-12+1):end);
    
%     sp_mean=mean(Clim_Sml_p.bio(:,lyr),2);
%     sf_mean=mean(Clim_Sml_f.bio(:,lyr),2);
%     sd_mean=mean(Clim_Sml_d.bio(:,lyr),2);
%     mp_mean=mean(Clim_Med_p.bio(:,lyr),2);
%     mf_mean=mean(Clim_Med_f.bio(:,lyr),2);
%     md_mean=mean(Clim_Med_d.bio(:,lyr),2);
%     lp_mean=mean(Clim_Lrg_p.bio(:,lyr),2);
%     ld_mean=mean(Clim_Lrg_d.bio(:,lyr),2);
%     b_mean =mean(Clim_Bent.bio(:,lyr),2);
    
    [ni,nt] = size(Clim_Lrg_d.yield);
    nyr = nt/12;
    mos = repmat(MNTH,ni,nyr);
    area_mat = repmat(area_km2,1,nt);
    
    MF.catch = Clim_Med_f.yield .*mos .*area_mat;
    MP.catch = Clim_Med_p.yield .*mos .*area_mat;
    MD.catch = Clim_Med_d.yield .*mos .*area_mat;
    LP.catch = Clim_Lrg_p.yield .*mos .*area_mat;
    LD.catch = Clim_Lrg_d.yield .*mos .*area_mat;
    
    % Total annual catch per grid cell
    mf_tc(:,n)=nansum(MF.catch(:,lyr),2);
    mp_tc(:,n)=nansum(MP.catch(:,lyr),2);
    md_tc(:,n)=nansum(MD.catch(:,lyr),2);
    lp_tc(:,n)=nansum(LP.catch(:,lyr),2);
    ld_tc(:,n)=nansum(LD.catch(:,lyr),2);
    
    % Total global annual catch
    mf_gtc(n)=nansum(MF.catch(:));
    mp_gtc(n)=nansum(MP.catch(:));
    md_gtc(n)=nansum(MD.catch(:));
    lp_gtc(n)=nansum(LP.catch(:));
    ld_gtc(n)=nansum(LD.catch(:));
    
end

save([fpath 'Climatol_fish_diffF_tot_catch.mat'],'Fnames','Fpaths','Files',...
    'mf_tc','mp_tc','md_tc','lp_tc','ld_tc',...
    'mf_gtc','mp_gtc','md_gtc','lp_gtc','ld_gtc');
