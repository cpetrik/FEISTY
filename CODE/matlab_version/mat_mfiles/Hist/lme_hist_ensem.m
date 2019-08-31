function [lme_mcatch,lme_mbio,lme_area] = lme_hist_ensem(sf_mean,sp_mean,sd_mean,...
    mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean,mf_my,mp_my,md_my,lp_my,ld_my,...
    fname,simname)


% Calc LME biomass of FEISTY
% Forecast time period 2006-2100
% 95 years
% Saved as mat files

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t','AREA_OCN');
grid = csvread([cpath 'grid_csv.csv']);
load([cpath 'lme_mask_esm2m.mat']);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);
tlme = lme_mask_esm2m';

%% Plots in space
[ni,nj]=size(geolon_t);

Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);

Bsf = Zmf;
Bsp = Zmf;
Bsd = Zmf;
Bmf = Zmf;
Bmp = Zmf;
Bmd = Zmf;
Blp = Zmf;
Bld = Zmf;
Bb = Zmf;

Zmf(grid(:,1))=mf_my;
Zmp(grid(:,1))=mp_my;
Zmd(grid(:,1))=md_my;
Zlp(grid(:,1))=lp_my;
Zld(grid(:,1))=ld_my;

Bsf(grid(:,1))=sf_mean;
Bsp(grid(:,1))=sp_mean;
Bsd(grid(:,1))=sd_mean;
Bmf(grid(:,1))=mf_mean;
Bmp(grid(:,1))=mp_mean;
Bmd(grid(:,1))=md_mean;
Blp(grid(:,1))=lp_mean;
Bld(grid(:,1))=ld_mean;
Bb(grid(:,1))=b_mean;

% g/m2 --> total g
Amf_mcatch = Zmf .* AREA_OCN * 365; %mean fish catch per yr
Amp_mcatch = Zmp .* AREA_OCN * 365;
Amd_mcatch = Zmd .* AREA_OCN * 365;
Alp_mcatch = Zlp .* AREA_OCN * 365;
Ald_mcatch= Zld .* AREA_OCN * 365;

Asf_mean = Bsf .* AREA_OCN;
Asp_mean = Bsp .* AREA_OCN;
Asd_mean = Bsd .* AREA_OCN;
Amf_mean = Bmf .* AREA_OCN;
Amp_mean = Bmp .* AREA_OCN;
Amd_mean = Bmd .* AREA_OCN;
Alp_mean = Blp .* AREA_OCN;
Ald_mean = Bld .* AREA_OCN;
Ab_mean = Bb .* AREA_OCN;

%% Calc LMEs
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
save([fname '_LME_' simname '.mat'],...
    'lme_mcatch','lme_mbio','lme_area');

end