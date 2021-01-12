% Plot effective TEs at LME scale
% Climatology
% 150 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
cdir='/Volumes/FEISTY/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);

%% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
% ENTER -100 TO MAP ORIGIN LONG

AREA_OCN = max(area,1);

%% FEISTY
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A080_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
dpath=['/Volumes/FEISTY/NC/Clim_comp_tests/' cfile '/NoNuUpdate_'];
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/clim_complete/post_proc/pp_figs/',...
    cfile,'/NoNuUpdate_'];

load([dpath 'TEeffDet_Climatol_All_fish03_' cfile '.mat']);

%% Calc LMEs
tlme = lme_mask_onedeg;

lme_te = NaN*ones(66,4);
for L=1:66
    lid = find(tlme==L);
    %TEeff
    lme_te(L,1) = nan;
    lme_te(L,2) = nanmean(TEeff_ATL(lid));
    lme_te(L,3) = nanmean(TEeff_HTL(lid));
    lme_te(L,4) = nanmean(TEeff_LTL(lid));
    
end

lme_m = NaN*ones(ni,nj);
lme_atl = lme_m;
lme_htl = lme_m;
lme_ltl = lme_m;

for L=1:66
    lid = find(tlme==L);

    lme_m(lid)      = lme_te(L,1);
    lme_atl(lid)    = lme_te(L,2);
    lme_htl(lid)    = lme_te(L,3);
    lme_ltl(lid)    = lme_te(L,4);
end

%%
save([dpath 'TEeffDet_Climatol_All_fish03_' cfile '.mat'],'lme_te',...
    'lme_atl','lme_htl','lme_ltl','-append');

TEATL = real(lme_te(:,2).^(1/4));
TEHTL = real(lme_te(:,3).^(1/3));
TELTL = real(lme_te(:,4).^(1/1.33333));

Tab=table([1:66]',lme_te(:,2),lme_te(:,3),lme_te(:,4),...
    TEATL,TEHTL,TELTL,...
    'VariableNames',{'LME','TEeffATL','TEeffHTL','TEeffLTL','TEATL',...
    'TEHTL','TELTL'});
writetable(Tab,[dpath 'LME_TEeff_clim_fished_',harv,'_' cfile '.csv'],'Delimiter',',');
save([dpath 'LME_TEeff_clim_fished_',harv,'_' cfile '.mat'],'Tab');


