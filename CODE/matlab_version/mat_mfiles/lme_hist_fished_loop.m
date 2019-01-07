% Calc LME biomass of POEM
% Hindcast time period 1861-2005
% 145 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t','AREA_OCN');
grid = csvread([cpath 'grid_csv.csv']);
load([cpath 'lme_mask_esm2m.mat']);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

kays = 0.0405:0.01:0.0916;
bees = 0.125:0.025:0.25; %bees = 0.1:0.05:0.35;
harv = '03';

for g = 1:length(bees)
    bpow = bees(g);
    
    for n = 1:length(kays)
        kt = kays(n);
        
        tkfn = num2str(100+int64(100*kt));
        tbfn = num2str(100+int64(100*bpow));
        
        cfile = ['Dc_enc70_cmax-metab20_b',tbfn(2:end),'_k',tkfn(2:end),'_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100'];
        
        fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
        ppath = [pp cfile '/'];
        if (~isdir(ppath))
            mkdir(ppath)
        end
        
        load([fpath 'Means_Hindcast_fished_' harv '_' cfile '.mat']);
        
        
        %% Plots in space
        [ni,nj]=size(geolon_t);
        
        Zmf=NaN*ones(ni,nj);
        Zmp=NaN*ones(ni,nj);
        Zmd=NaN*ones(ni,nj);
        Zlp=NaN*ones(ni,nj);
        Zld=NaN*ones(ni,nj);
        
        Zmf(grid(:,1))=mf_my;
        Zmp(grid(:,1))=mp_my;
        Zmd(grid(:,1))=md_my;
        Zlp(grid(:,1))=lp_my;
        Zld(grid(:,1))=ld_my;
        
        % g/m2 --> total g
        Amf_mcatch = Zmf .* AREA_OCN * 365; %mean fish catch per yr
        Amp_mcatch = Zmp .* AREA_OCN * 365;
        Amd_mcatch = Zmd .* AREA_OCN * 365;
        Alp_mcatch = Zlp .* AREA_OCN * 365;
        Ald_mcatch= Zld .* AREA_OCN * 365;
        
        
        %% Calc LMEs
        lat = geolat_t;
        lon = geolon_t;
        
        tlme = lme_mask_esm2m';
        
        lme_mcatch = NaN*ones(66,5);
        
        for L=1:66
            lid = find(tlme==L);
            %total catch g
            lme_mcatch(L,1) = nansum(Amf_mcatch(lid));
            lme_mcatch(L,2) = nansum(Amp_mcatch(lid));
            lme_mcatch(L,3) = nansum(Amd_mcatch(lid));
            lme_mcatch(L,4) = nansum(Alp_mcatch(lid));
            lme_mcatch(L,5) = nansum(Ald_mcatch(lid));
            
        end
        
        %%
        save([fpath 'LME_hindcast_fished',harv,'_' cfile '.mat'],...
            'lme_mcatch');
        
    end
end