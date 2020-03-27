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

%load([dpath 'TEeffDetZprod_Historic_All_fish03_' cfile '.mat']);

%% Calc TEs
load([dpath 'Means_Historic_' harv '_prod_' cfile '.mat']);

geolon_t=double(geolon_t);
geolat_t=double(geolat_t);
[ni,nj]=size(geolon_t);

Psf=NaN*ones(ni,nj);
Psp=NaN*ones(ni,nj);
Psd=NaN*ones(ni,nj);
Pmf=NaN*ones(ni,nj);
Pmp=NaN*ones(ni,nj);
Pmd=NaN*ones(ni,nj);
Plp=NaN*ones(ni,nj);
Pld=NaN*ones(ni,nj);
Pb=NaN*ones(ni,nj);

Psf(ID)=sf_prod50;
Psp(ID)=sp_prod50;
Psd(ID)=sd_prod50;
Pmf(ID)=mf_prod50;
Pmp(ID)=mp_prod50;
Pmd(ID)=md_prod50;
Plp(ID)=lp_prod50;
Pld(ID)=ld_prod50;
Pb(ID)=b_mean50;

Psf(Psf(:)<0) = 0;
Psp(Psp(:)<0) = 0;
Psd(Psd(:)<0) = 0;
Pmf(Pmf(:)<0) = 0;
Pmp(Pmp(:)<0) = 0;
Pmd(Pmd(:)<0) = 0;
Plp(Plp(:)<0) = 0;
Pld(Pld(:)<0) = 0;
Pb(Pb(:)<0) = 0;

All = Psp+Psf+Psd+Pmp+Pmf+Pmd+Plp+Pld;
AllF = Psf+Pmf;
AllP = Psp+Pmp+Plp;
AllD = Psd+Pmd+Pld;
AllS = Psp+Psf+Psd;
AllM = Pmp+Pmf+Pmd;
AllL = Plp+Pld;

%% Zoop and det and npp 

%ESM2M in mmol N m-2 or mmol N m-2 d-1
% molN/m2 --> g/m2
% 106/16 mol C in 1 mol N
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W
mz_mean_hist = mz_mean_hist * (106.0/16.0) * 12.01 * 9.0;
lz_mean_hist = lz_mean_hist * (106.0/16.0) * 12.01 * 9.0;
% molN/m2/s --> g/m2/d
mzprod_mean_hist = mzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_mean_hist = lzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_mean_hist = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_hist = npp_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

% molN/m2/s --> gC
mnpp = npp_mean_hist * (106.0/16.0) * 12.01 * 60 * 60 * 24 * 365 .* AREA_OCN;

mmz_mean = mz_mean_hist;
mlz_mean = lz_mean_hist;
mmz_prod = mzprod_mean_hist;
mlz_prod = lzprod_mean_hist;
mdet = det_mean_hist;

%% Effective TEs
% With BE*det instead of Bent
BE = 0.075;

teS = AllS./(mmz_prod);
teS(teS==-Inf) = NaN;
teS(teS==Inf) = NaN;
teS(teS<0) = NaN;

teM = AllM./(BE*mdet + mmz_prod + mlz_prod + AllS);
teM(teM==-Inf) = NaN;
teM(teM==Inf) = NaN;
teM(teM<0) = NaN;

teL = AllL./(BE*mdet + AllM);
teL(teL==-Inf) = NaN;
teL(teL==Inf) = NaN;
teL(teL<0) = NaN;

%TEeff_ATL = production_L/NPP
TEeff_ATL = AllL ./ npp_hist;
TEeff_ATL(TEeff_ATL==-Inf) = NaN;
TEeff_ATL(TEeff_ATL==Inf) = NaN;
TEeff_ATL(TEeff_ATL<0) = NaN;

%TEeff_LTL = (production_benthic_invert+mesozoo_prod_to_fish)/NPP
TEeff_LTL = (BE*mdet + mmz_prod + mlz_prod) ./ npp_hist;
TEeff_LTL(TEeff_LTL==-Inf) = NaN;
TEeff_LTL(TEeff_LTL==Inf) = NaN;
TEeff_LTL(TEeff_LTL<0) = NaN;

%TEeff_HTL = production_L/(production_benthic_invert+mesozoo_prod_to_fish)
TEeff_HTL = AllL ./ (BE*mdet + mmz_prod + mlz_prod); 
TEeff_HTL(TEeff_HTL==-Inf) = NaN;
TEeff_HTL(TEeff_HTL==Inf) = NaN;
TEeff_HTL(TEeff_HTL<0) = NaN;

% TELTL = real(TEeff_LTL.^(1/1.333)); %(1+1+2/3)
% TEATL = real(TEeff_ATL.^(1/4));     %should this be 1/3?
% TEHTL = real(TEeff_HTL.^(1/3));     %should this be 1/2?

%% Calc LMEs
tlme = lme_mask_esm2m';

for L=1:66
    lid = find(tlme==L);
    %TEeff
    lme_teS(L,1) = nanmean(teS(lid));
    lme_teM(L,1) = nanmean(teM(lid));
    lme_teL(L,1) = nanmean(teL(lid));
    lme_teATL(L,1) = nanmean(TEeff_ATL(lid));
    lme_teHTL(L,1) = nanmean(TEeff_HTL(lid));
    lme_teLTL(L,1) = nanmean(TEeff_LTL(lid));
    lme_npp(L,1) = nansum(mnpp(lid));
end

lme_m = NaN*ones(ni,nj);
lme_l = lme_m;
lme_s = lme_m;
lme_htl = lme_m;
lme_ltl = lme_m;
lme_atl = lme_m;
for L=1:66
    lid = find(tlme==L);
    
    lme_s(lid)      = lme_teS(L);
    lme_m(lid)      = lme_teM(L);
    lme_l(lid)      = lme_teL(L);
    lme_atl(lid)   = lme_teATL(L);
    lme_htl(lid)   = lme_teHTL(L);
    lme_ltl(lid)   = lme_teLTL(L);
end

save([dpath 'TEms_LMEs_TEeffDetZprod_Historic_All_fish03_' cfile '.mat'],'lme_s',...
    'lme_m','lme_l','lme_htl','lme_ltl','lme_atl');

%%
alme = 1:66;
nlm1 = find(isnan(tlme)); %non-lme areas, 
nupw = find(abs(geolat_t(:))>10); %non-upwelling
nlme = intersect(nlm1,nupw); %MAYBE CHANGE THIS TO INCLUDE EQ
ulme = [3,13,27,29];
clme = setdiff(alme,ulme);

%MAYBE TAKE QUANTILE OF ALL CELLS IN LMES VS. QUANTILE OF LME MEANS

biome_mmte = NaN*ones(12,3);
biome_npp = NaN*ones(3,1);
for n=1:2
    if n==1
        L=clme;
    elseif n==2
        L=ulme;
    end
    biome_mmte(1,n) = quantile(lme_teS(L),0.01);
    biome_mmte(2,n) = quantile(lme_teS(L),0.99);
    
    biome_mmte(3,n) = quantile(lme_teM(L),0.01);
    biome_mmte(4,n) = quantile(lme_teM(L),0.99);
    
    biome_mmte(5,n) = quantile(lme_teL(L),0.01);
    biome_mmte(6,n) = quantile(lme_teL(L),0.99);
    
    biome_mmte(7,n) = quantile(lme_teLTL(L),0.01);
    biome_mmte(8,n) = quantile(lme_teLTL(L),0.99);
    
    biome_mmte(9,n) = quantile(lme_teHTL(L),0.01);
    biome_mmte(10,n) = quantile(lme_teHTL(L),0.99);
    
    biome_mmte(11,n) = quantile(lme_teATL(L),0.01);
    biome_mmte(12,n) = quantile(lme_teATL(L),0.99);
    
    biome_npp(n) = nansum(lme_npp(L));
end

biome_mmte(1,3) = quantile(teS(nlme),0.01);
biome_mmte(2,3) = quantile(teS(nlme),0.99);
biome_mmte(3,3) = quantile(teM(nlme),0.01);
biome_mmte(4,3) = quantile(teM(nlme),0.99);
biome_mmte(5,3) = quantile(teL(nlme),0.01);
biome_mmte(6,3) = quantile(teL(nlme),0.99);
biome_mmte(7,3) = quantile(TEeff_LTL(nlme),0.01);
biome_mmte(8,3) = quantile(TEeff_LTL(nlme),0.99);
biome_mmte(9,3) = quantile(TEeff_HTL(nlme),0.01);
biome_mmte(10,3) = quantile(TEeff_HTL(nlme),0.99);
biome_mmte(11,3) = quantile(TEeff_ATL(nlme),0.01);
biome_mmte(12,3) = quantile(TEeff_ATL(nlme),0.99);

biome_npp(3) = nansum(mnpp(nlm1));

%%
TELTL = real(lme_teLTL.^(1/1.333)); %(1+1+2/3)
TEATL = real(lme_teATL.^(1/4));     
TEHTL = real(lme_teHTL.^(1/3));

Tab=table([1:66]',lme_teS,lme_teM,lme_teL,TELTL,TEHTL,TEATL,...
    'VariableNames',{'LME','TES','TEM','TEL','TELTL','TEHTL','TEATL'});
writetable(Tab,[dpath 'TEms_LMEs_TEeffDetZprod_hist_',harv,'_' cfile '.csv'],'Delimiter',',','WriteRowNames',true);
save([dpath 'TEms_LMEs_TEeffDetZprod_Historic_All_fish03_' cfile '.mat'],...
    'Tab','-append');

%% Conversions to TE and save
biome_mm(1:6,:) = biome_mmte(1:6,:);
biome_mm(7:8,:) = real(biome_mmte(7:8,:).^(1/1.3333)); %TELTL
biome_mm(9:10,:) = real(biome_mmte(9:10,:).^(1/3));   %TEHTL should this be 1/2?
biome_mm(11:12,:) = real(biome_mmte(11:12,:).^(1/4));   %TEATL should this be 1/3?

Q = array2table(biome_mm,'VariableNames',{'Coastal','Upwelling','Oceanic'},...
    'RowNames',{'minTES','maxTES','minTEM','maxTEM','minTEL','maxTEL',...
    'minTELTL','maxTELTL','minTEHTL','maxTEHTL','minTEATL','maxTEATL'});

B = array2table(biome_mmte,'VariableNames',{'Coastal','Upwelling','Oceanic'},...
    'RowNames',{'minTES','maxTES','minTEM','maxTEM','minTEL','maxTEL',...
    'minTEeffLTL','maxTEeffLTL','minTEeffHTL','maxTEeffHTL','minTEeffATL',...
    'maxTEeffATL'});

writetable(Q,[dpath 'TEms_Ryther_TEDetZprod_hist_',harv,'_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(B,[dpath 'TEms_Ryther_TEeffDetZprod_hist_',harv,'_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
save([dpath 'TEms_Ryther_TEeffDetZprod_hist_',harv,'_' cfile '.mat'],'Q','B');


%% TAKE QUANTILE OF ALL CELLS IN LMES VS. QUANTILE OF LME MEANS

lid = cell(66,1);
for L=1:66
    lid{L,1} = find(tlme==L);
end

uid = [];
for n=1:length(ulme)
    L = ulme(n);
    uid = [uid; lid{L,1}];
end

cid = [];
for n=1:length(clme)
    L = clme(n);
    cid = [cid; lid{L,1}];
end

%%
biome_mTE = NaN*ones(12,3);
biome_mNPP = NaN*ones(3,1);
biome_teATL = NaN*ones(7,3);
for n=1:3
    if n==1
        L=cid;
        Ln=cid;
    elseif n==2
        L=uid;
        Ln=uid;
    elseif n==3
        L=nlme;
        Ln=nlm1;
    end
    biome_mTE(1,n) = quantile(teS(L),0.01);
    biome_mTE(2,n) = quantile(teS(L),0.99);
    biome_mTE(3,n) = quantile(teM(L),0.01);
    biome_mTE(4,n) = quantile(teM(L),0.99);
    biome_mTE(5,n) = quantile(teL(L),0.01);
    biome_mTE(6,n) = quantile(teL(L),0.99);
    biome_mTE(7,n) = quantile(TEeff_LTL(L),0.01);
    biome_mTE(8,n) = quantile(TEeff_LTL(L),0.99);
    biome_mTE(9,n) = quantile(TEeff_HTL(L),0.01);
    biome_mTE(10,n) = quantile(TEeff_HTL(L),0.99);
    biome_mTE(11,n) = quantile(TEeff_ATL(L),0.01);
    biome_mTE(12,n) = quantile(TEeff_ATL(L),0.99);
    
    biome_teATL(1,n) = quantile(TEeff_ATL(L),0.01);
    biome_teATL(2,n) = quantile(TEeff_ATL(L),0.05);
    biome_teATL(3,n) = quantile(TEeff_ATL(L),0.10);
    biome_teATL(4,n) = quantile(TEeff_ATL(L),0.50);
    biome_teATL(5,n) = quantile(TEeff_ATL(L),0.90);
    biome_teATL(6,n) = quantile(TEeff_ATL(L),0.95);
    biome_teATL(7,n) = quantile(TEeff_ATL(L),0.99);
    
%     biome_mTE(1,n) = quantile(teS(L),0.1);
%     biome_mTE(2,n) = quantile(teS(L),0.9);
%     biome_mTE(3,n) = quantile(teM(L),0.1);
%     biome_mTE(4,n) = quantile(teM(L),0.9);
%     biome_mTE(5,n) = quantile(teL(L),0.1);
%     biome_mTE(6,n) = quantile(teL(L),0.9);
%     biome_mTE(7,n) = quantile(TEeff_LTL(L),0.1);
%     biome_mTE(8,n) = quantile(TEeff_LTL(L),0.9);
%     biome_mTE(9,n) = quantile(TEeff_HTL(L),0.1);
%     biome_mTE(10,n) = quantile(TEeff_HTL(L),0.9);
%     biome_mTE(11,n) = quantile(TEeff_ATL(L),0.1);
%     biome_mTE(12,n) = quantile(TEeff_ATL(L),0.9);
    
%     biome_mTE(1,n) = quantile(teS(L),0.05);
%     biome_mTE(2,n) = quantile(teS(L),0.95);
%     biome_mTE(3,n) = quantile(teM(L),0.05);
%     biome_mTE(4,n) = quantile(teM(L),0.95);
%     biome_mTE(5,n) = quantile(teL(L),0.05);
%     biome_mTE(6,n) = quantile(teL(L),0.95);
%     biome_mTE(7,n) = quantile(TEeff_LTL(L),0.05);
%     biome_mTE(8,n) = quantile(TEeff_LTL(L),0.95);
%     biome_mTE(9,n) = quantile(TEeff_HTL(L),0.05);
%     biome_mTE(10,n) = quantile(TEeff_HTL(L),0.95);
%     biome_mTE(11,n) = quantile(TEeff_ATL(L),0.05);
%     biome_mTE(12,n) = quantile(TEeff_ATL(L),0.95);

    biome_mNPP(n) = nansum(mnpp(Ln));
end

%% Conversions to TE and save
mbiome(1:6,:) = biome_mTE(1:6,:);
mbiome(7:8,:) = real(biome_mTE(7:8,:).^(1/1.3333)); %TELTL
mbiome(9:10,:) = real(biome_mTE(9:10,:).^(1/3));    %TEHTL should this be 1/2?
mbiome(11:12,:) = real(biome_mTE(11:12,:).^(1/4));  %TEATL should this be 1/3?

biome_mATL = real(biome_teATL.^(1/4));  %TEATL should this be 1/3?

Q2 = array2table(mbiome,'VariableNames',{'Coastal','Upwelling','Oceanic'},...
    'RowNames',{'minTES','maxTES','minTEM','maxTEM','minTEL','maxTEL',...
    'minTELTL','maxTELTL','minTEHTL','maxTEHTL','minTEATL','maxTEATL'});

B2 = array2table(biome_mTE,'VariableNames',{'Coastal','Upwelling','Oceanic'},...
    'RowNames',{'minTES','maxTES','minTEM','maxTEM','minTEL','maxTEL',...
    'minTEeffLTL','maxTEeffLTL','minTEeffHTL','maxTEeffHTL','minTEeffATL',...
    'maxTEeffATL'});

A2 = array2table(biome_mATL,'VariableNames',{'Coastal','Upwelling','Oceanic'},...
    'RowNames',{'1','5','10','50','90','95','99'});

writetable(Q2,[dpath 'TEms_Ryther_TEDetZprod_1-99_hist_',harv,'_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(B2,[dpath 'TEms_Ryther_TEeffDetZprod_1-99_hist_',harv,'_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(A2,[dpath 'TEms_Ryther_TEeffDetZprod_ATLquant_hist_',harv,'_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
save([dpath 'TEms_Ryther_TEeffDetZprod_1-99_hist_',harv,'_' cfile '.mat'],'Q2','B2');


NPP = array2table(biome_mNPP,'RowNames',{'Coastal','Upwelling','Oceanic'});
writetable(NPP,[dpath 'TEms_Ryther_NPP_hist_',harv,'_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
