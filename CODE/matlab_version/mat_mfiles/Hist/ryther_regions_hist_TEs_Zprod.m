% Plot effective TEs at LME scale
% Historic
% 1860-2005, last 50 years
% Saved as mat files
% Use Zprod instead of loss

clear all
close all

dp = '/Volumes/FEISTY/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
Pdir = '/Volumes/FEISTY/POEM_JLD/esm2m_hist/';
cdir = '/Volumes/FEISTY/GCM_DATA/ESM2M_hist/';
load([gpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t','AREA_OCN');
load([gpath 'lme_mask_esm2m.mat']);
grid = csvread([gpath 'grid_csv.csv']);
load([cpath 'cobalt_det_temp_zoop_npp_means.mat']);
%lme areas?

geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
[ni,nj]=size(geolon_t);
AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

ID = grid(:,1);

%%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'TEeffDetZprod_Historic_All_fish03_' cfile '.mat']);

%% Calc LMEs
tlme = lme_mask_esm2m';

lme_te = NaN*ones(66,2);
for L=1:66
    lid = find(tlme==L);
    %TEeff
    lme_te(L,1) = nanmean(TEeffM(lid));
    lme_te(L,2) = nanmean(TEeff_L(lid));
    lme_te(L,4) = nanmean(TEeff_HTLd(lid));
    lme_te(L,6) = nanmean(TEeff_LTLd(lid));
    
end

lme_m = NaN*ones(ni,nj);
lme_l = lme_m;
lme_htlD = lme_m;
lme_ltlD = lme_m;
for L=1:66
    lid = find(tlme==L);
    
    lme_m(lid)      = lme_te(L,1);
    lme_l(lid)      = lme_te(L,2);
    lme_htlD(lid)   = lme_te(L,4);
    lme_ltlD(lid)   = lme_te(L,6);
end

save([dpath 'TEeffDetZprod_Historic_All_fish03_' cfile '.mat'],'lme_te',...
    'lme_m','lme_l','lme_htlD','lme_ltlD','-append');

%%
alme = 1:66;
nlm1 = find(isnan(tlme)); %non-lme areas, 
nupw = find(abs(geolat_t(:))>10); %non-upwelling
nlme = intersect(nlm1,nupw);
ulme = [3,13,27,29];
clme = setdiff(alme,ulme);

biome_mmte = NaN*ones(8,3);
for n=1:2
    if n==1
        L=clme;
    elseif n==2
        L=ulme;
    end
    biome_mmte(1,n) = quantile(lme_te(L,1),0.01);
    biome_mmte(2,n) = quantile(lme_te(L,1),0.99);
    biome_mmte(3,n) = quantile(lme_te(L,2),0.01);
    biome_mmte(4,n) = quantile(lme_te(L,2),0.99);
    biome_mmte(5,n) = quantile(lme_te(L,4),0.01);
    biome_mmte(6,n) = quantile(lme_te(L,4),0.99);
    biome_mmte(7,n) = quantile(lme_te(L,6),0.01);
    biome_mmte(8,n) = quantile(lme_te(L,6),0.99);
end

biome_mmte(1,3) = quantile(TEeffM(nlme),0.01);
biome_mmte(2,3) = quantile(TEeffM(nlme),0.99);
biome_mmte(3,3) = quantile(TEeff_L(nlme),0.01);
biome_mmte(4,3) = quantile(TEeff_L(nlme),0.99);
biome_mmte(5,3) = quantile(TEeff_HTLd(nlme),0.01);
biome_mmte(6,3) = quantile(TEeff_HTLd(nlme),0.99);
biome_mmte(7,3) = quantile(TEeff_LTLd(nlme),0.01);
biome_mmte(8,3) = quantile(TEeff_LTLd(nlme),0.99);

%%
TEM   = real(lme_te(:,1).^(1/2));
TEATL = real(lme_te(:,2).^(1/4));
TEHTL = real(lme_te(:,4).^(1/3));
TELTL = real(lme_te(:,6).^(1/1.3333));

Tab=table([1:66]',lme_te(:,2),lme_te(:,4),lme_te(:,6),...
    TEATL,TEHTL,TELTL,...
    'VariableNames',{'LME','TEeffATL','TEeffHTL','TEeffLTL',...
    'TEATL','TEHTL','TELTL'});
writetable(Tab,[dpath 'LME_TEeffDetZprod_hist_',harv,'_' cfile '.csv'],'Delimiter',',','WriteRowNames',true);
save([dpath 'LME_TEeff_hist_',harv,'_' cfile '.mat'],'Tab');

% Conversions to TE and save
biome_mm(1:2,:) = real(biome_mmte(1:2,:).^(1/2));   %TEM should this be 1/1?
biome_mm(3:4,:) = real(biome_mmte(3:4,:).^(1/4));   %TEL should this be 1/3?
biome_mm(5:6,:) = real(biome_mmte(5:6,:).^(1/3));   %TEHTL should this be 1/2?
biome_mm(7:8,:) = real(biome_mmte(7:8,:).^(1/1.3333));

Q = array2table(biome_mm,'VariableNames',{'Coastal','Upwelling','Oceanic'},...
    'RowNames',{'minTEM','maxTEM','minTEATL','maxTEATL','minTEHTL','maxTEHTL',...
    'minTELTL','maxTELTL'});

B = array2table(biome_mmte,'VariableNames',{'Coastal','Upwelling','Oceanic'},...
    'RowNames',{'minTEeffM','maxTEeffM','minTEeffATL','maxTEeffATL','minTEeffHTL',...
    'maxTEeffHTL','minTEeffLTL','maxTEeffLTL'});

writetable(Q,[dpath 'Ryther_TEDetZprod_hist_',harv,'_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(B,[dpath 'Ryther_TEeffDetZprod_hist_',harv,'_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
save([dpath 'Ryther_TEeffDetZprod_hist_',harv,'_' cfile '.mat'],'Q','B');


%% repeat with 90-95 means
load([dpath 'TEeffDet_Historic9095_All_fish03_' cfile '.mat']);
load([cpath 'cobalt_hist9095_det_temp_zoop_npp_means.mat']);