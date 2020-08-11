function [lme_mcatch,lme_area] = lme_hist_fished_1meso_loop(mf_my50,...
    mp_my50,md_my50,lp_my50,ld_my50,dpath)

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';

load([gpath 'lme_mask_esm2m.mat'],'lme_mask_esm2m');
load([gpath 'hindcast_gridspec.mat'],'AREA_OCN','geolon_t');
grid = csvread([gpath 'grid_csv.csv']);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);
tlme = lme_mask_esm2m';

harv = 'All_fish03';

%% Plots in space
[ni,nj]=size(geolon_t);

Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);

Zmf(grid(:,1))=mf_my50;
Zmp(grid(:,1))=mp_my50;
Zmd(grid(:,1))=md_my50;
Zlp(grid(:,1))=lp_my50;
Zld(grid(:,1))=ld_my50;

% g/m2 --> total g
Amf_mcatch = Zmf .* AREA_OCN * 365; %mean fish catch per yr
Amp_mcatch = Zmp .* AREA_OCN * 365;
Amd_mcatch = Zmd .* AREA_OCN * 365;
Alp_mcatch = Zlp .* AREA_OCN * 365;
Ald_mcatch= Zld .* AREA_OCN * 365;

%% Calc LMEs
lme_mcatch = NaN*ones(66,5);
lme_area = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    %total catch g
    lme_mcatch(L,1) = nansum(Amf_mcatch(lid));
    lme_mcatch(L,2) = nansum(Amp_mcatch(lid));
    lme_mcatch(L,3) = nansum(Amd_mcatch(lid));
    lme_mcatch(L,4) = nansum(Alp_mcatch(lid));
    lme_mcatch(L,5) = nansum(Ald_mcatch(lid));
    
    %total area of LME
    lme_area(L,1) = nansum(AREA_OCN(lid));
end

%%
save([dpath 'LME_Hist_1meso_',harv,'_.mat'],...
    'lme_mcatch','lme_area');



