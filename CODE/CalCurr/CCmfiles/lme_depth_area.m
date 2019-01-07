% Calc frac of LME area <200 m 
% Climatology

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);

%%
AREA_OCN = max(area,1);
[ni,nj]=size(lon);
tlme = lme_mask_onedeg;
%%
lme_area = zeros(66,2);
lme_depth = zeros(66,1);
for L=1:66
    lid = find(tlme==L);
    dep = depth(lid);
    shal = (dep<200);
    ar = AREA_OCN(lid);
    %total area <200 and greater/equal
    if (~isempty(ar(shal==1)))
        lme_area(L,1) = nansum(ar(shal==1));
    end
    %area-weighted mean depth
    if (~isempty(ar(shal==0)))
        lme_area(L,2) = nansum(ar(shal==0));
    end
    lme_depth(L,1) = nansum(dep.*ar)./nansum(ar);
end

lme_shal_frac = lme_area(:,1)./(lme_area(:,1)+lme_area(:,2));

%%


