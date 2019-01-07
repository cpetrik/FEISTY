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

AREA_OCN = max(area,1);
tlme = lme_mask_onedeg;

%%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

ktemp = 0.0805:0.005:0.1005;
bees = 0.025:0.025:0.125;

lme_mcatch_mf = NaN*ones(length(bees),length(ktemp),66);
lme_mcatch_mp = lme_mcatch_mf;
lme_mcatch_md = lme_mcatch_mf;
lme_mcatch_lp = lme_mcatch_mf;
lme_mcatch_ld = lme_mcatch_mf;
lme_mbio_sf = NaN*ones(length(bees),length(ktemp),66);
lme_mbio_sp = lme_mbio_sf;
lme_mbio_sd = lme_mbio_sf;
lme_mbio_mf = lme_mbio_sf;
lme_mbio_mp = lme_mbio_sf;
lme_mbio_md = lme_mbio_sf;
lme_mbio_lp = lme_mbio_sf;
lme_mbio_ld = lme_mbio_sf;
lme_mbio_b = lme_mbio_sf;
lme_area = NaN*ones(66,1);

for k=1:length(ktemp)
    for j=1:length(bees)   
        bent_eff = bees(j); 
        kt = ktemp(k);   
        tkfn = num2str(1000+int64(1000*kt));
        tbe = num2str(100+int64(100*bent_eff));
        
        sfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,...
            '/param_sens/Climatol_All_fish03__m-k',tkfn(2:end),'_BE',tbe(2:end),'.mat'];
           %'/Climatol_All_fish03__m-k096_BE10.mat'];
        load(sfile);
        
        %% Plots in space
        [ni,nj]=size(lon);
        
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
        
        Zsf(ID)=sf_mean;
        Zsp(ID)=sp_mean;
        Zsd(ID)=sd_mean;
        Zmf(ID)=mf_mean;
        Zmp(ID)=mp_mean;
        Zmd(ID)=md_mean;
        Zlp(ID)=lp_mean;
        Zld(ID)=ld_mean;
        Zb(ID)=b_mean;
        
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
        % g/m2 --> total g
        Asf_mean = Zsf .* AREA_OCN;
        Asp_mean = Zsp .* AREA_OCN;
        Asd_mean = Zsd .* AREA_OCN;
        Amf_mean = Zmf .* AREA_OCN;
        Amp_mean = Zmp .* AREA_OCN;
        Amd_mean = Zmd .* AREA_OCN;
        Alp_mean = Zlp .* AREA_OCN;
        Ald_mean = Zld .* AREA_OCN;
        Ab_mean  = Zb .* AREA_OCN;
        
        %% Calc LMEs
        
        for L=1:66
            lid = find(tlme==L);
            %total catch g
            lme_mcatch_mf(j,k,L) = nansum(Amf_mcatch(lid));
            lme_mcatch_mp(j,k,L) = nansum(Amp_mcatch(lid));
            lme_mcatch_md(j,k,L) = nansum(Amd_mcatch(lid));
            lme_mcatch_lp(j,k,L) = nansum(Alp_mcatch(lid));
            lme_mcatch_ld(j,k,L) = nansum(Ald_mcatch(lid));
            %mean biomass
            lme_mbio_sf(j,k,L) = nanmean(Asf_mean(lid));
            lme_mbio_sp(j,k,L) = nanmean(Asp_mean(lid));
            lme_mbio_sd(j,k,L) = nanmean(Asd_mean(lid));
            lme_mbio_mf(j,k,L) = nanmean(Amf_mean(lid));
            lme_mbio_mp(j,k,L) = nanmean(Amp_mean(lid));
            lme_mbio_md(j,k,L) = nanmean(Amd_mean(lid));
            lme_mbio_lp(j,k,L) = nanmean(Alp_mean(lid));
            lme_mbio_ld(j,k,L) = nanmean(Ald_mean(lid));
            lme_mbio_b(j,k,L) = nanmean(Ab_mean(lid));
            %total area of LME
            if (j==1 && k==1)
                lme_area(L,1) = nansum(AREA_OCN(lid));
            end
        end
        
        
    end
end
%%
save(['/Volumes/GFDL/NC/Matlab_new_size/',cfile,...
    '/param_sens/LME_Climatol_All_fish03_kt_BE_search.mat']);

