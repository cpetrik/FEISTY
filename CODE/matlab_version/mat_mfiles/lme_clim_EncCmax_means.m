% Calc LME biomass of POEM
% Climatology
% 100 years
% Only last 12 months saved to compare to SAUP for tuning

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/CSV/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cdir 'temp_100_1deg_ESM26_5yr_clim_191_195.mat'])
load([cdir 'btm_temp_1deg_ESM26_5yr_clim_191_195.mat'])

ptemp_mean_clim=squeeze(nanmean(temp_100,1));
btemp_mean_clim=squeeze(nanmean(btm_temp,1));

AREA_OCN = max(area,1);

encs = 10:10:100;
cmaxs = 10:10:50;
frate = 0.3;

%%
for c=4:length(cmaxs)
    for e=1:length(encs)
            gam = encs(e);
            h = cmaxs(c);
            tcfn = num2str(h);
            tefn = num2str(round(gam));
            
            cfile = ['Dc_enc',tefn,'-b200_cm',tcfn,'_m-b175-k09_fcrit20',...
                '_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100'];
            fpath=['/Volumes/GFDL/CSV/Matlab_new_size/' cfile '/'];
            
            ppath = [pp cfile '/'];
            dpath = [dp cfile '/'];
            
            harv = 'All_fish03';
            sfile = ['/Clim_means_All_fish03.mat'];
            load([dpath sfile]);
            
            %% Put into grid
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
            tlme = lme_mask_onedeg;
            
            lme_mcatch = NaN*ones(66,5);
            lme_mbio = NaN*ones(66,9);
            lme_area = NaN*ones(66,1);
            
            for L=1:66
                lid = find(tlme==L);
                %total catch g
                lme_mcatch(L,1) = nansum(Amf_mcatch(lid));
                lme_mcatch(L,2) = nansum(Amp_mcatch(lid));
                lme_mcatch(L,3) = nansum(Amd_mcatch(lid));
                lme_mcatch(L,4) = nansum(Alp_mcatch(lid));
                lme_mcatch(L,5) = nansum(Ald_mcatch(lid));
                %mean biomass
                lme_mbio(L,1) = nanmean(Asf_mean(lid));
                lme_mbio(L,2) = nanmean(Asp_mean(lid));
                lme_mbio(L,3) = nanmean(Asd_mean(lid));
                lme_mbio(L,4) = nanmean(Amf_mean(lid));
                lme_mbio(L,5) = nanmean(Amp_mean(lid));
                lme_mbio(L,6) = nanmean(Amd_mean(lid));
                lme_mbio(L,7) = nanmean(Alp_mean(lid));
                lme_mbio(L,8) = nanmean(Ald_mean(lid));
                lme_mbio(L,9) = nanmean(Ab_mean(lid));
                %total area of LME
                lme_area(L,1) = nansum(AREA_OCN(lid));
            end
            
            %%
            save([dpath 'LME_clim_',harv,'_loop.mat'],...
                'lme_mcatch','lme_mbio','lme_area');

    end
end

