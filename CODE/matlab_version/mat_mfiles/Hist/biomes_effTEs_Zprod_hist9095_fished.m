% Visualize output of FEISTY Historic globally
% 1990-1995, monthly means saved
% Transfer efficiency ("effective") 
% Use BE*det
% Fix TEeff estimates to use zoo prod instead of loss!!!

clear all
close all

Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
gpath='/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

% POEM
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
BE = 0.075;
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'TEeffDetZprod_Historic9095_All_fish03_' cfile '.mat']);

cmYOR=cbrewer('seq','YlOrRd',50,'PCHIP');
cmRP=cbrewer('seq','RdPu',50,'PCHIP');
cmPR=cbrewer('seq','PuRd',50,'PCHIP');

%% plot info

geolon_t=double(geolon_t);
geolat_t=double(geolat_t);
[ni,nj]=size(geolon_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

%% Effective TEs
% With BE*det instead of Bent

% Take means over biomes
load([gpath 'COBALT_biomes_ESM2M_9095.mat'])
tlme = biome_hist;

biome_te = NaN*ones(3,4);
biome_teeff = NaN*ones(3,4);
biome_mmte = NaN*ones(3,8);
for L=1:3
    lid = find(tlme==L);
    %TEeff
    biome_teeff(L,1) = nanmean(TEeffM(lid));
    biome_teeff(L,2) = nanmean(TEeff_L(lid));
    biome_teeff(L,3) = nanmean(TEeff_HTLd(lid));
    biome_teeff(L,4) = nanmean(TEeff_LTLd(lid));
    
    biome_te(L,1) = nanmean(real(TEeffM(lid).^(1/2)));
    biome_te(L,2) = nanmean(real(TEeff_L(lid).^(1/4)));
    biome_te(L,3) = nanmean(real(TEeff_HTLd(lid).^(1/3)));
    biome_te(L,4) = nanmean(real(TEeff_LTLd(lid).^(1/1.3333)));
    
    biome_mmte(L,1) = quantile(TEeffM(lid),0.01);
    biome_mmte(L,2) = quantile(TEeffM(lid),0.99);
    biome_mmte(L,3) = quantile(TEeff_L(lid),0.01);
    biome_mmte(L,4) = quantile(TEeff_L(lid),0.99);
    biome_mmte(L,5) = quantile(TEeff_HTLd(lid),0.01);
    biome_mmte(L,6) = quantile(TEeff_HTLd(lid),0.99);
    biome_mmte(L,7) = quantile(TEeff_LTLd(lid),0.01);
    biome_mmte(L,8) = quantile(TEeff_LTLd(lid),0.99);
end

biome_m = NaN*ones(ni,nj);
biome_l = biome_m;
biome_htlD = biome_m;
biome_ltlD = biome_m;
for L=1:3
    lid = find(tlme==L);

    biome_m(lid)      = biome_teeff(L,1);
    biome_l(lid)      = biome_teeff(L,2);
    biome_htlD(lid)   = biome_teeff(L,3);
    biome_ltlD(lid)   = biome_teeff(L,4);
end

%% Conversions to TE and save
biome_mm(:,1:2) = real(biome_mmte(:,1:2).^(1/2));   %TEM should this be 1/1?
biome_mm(:,3:4) = real(biome_mmte(:,3:4).^(1/4));   %TEL should this be 1/3?
biome_mm(:,5:6) = real(biome_mmte(:,5:6).^(1/3));   %TEHTL should this be 1/2?
biome_mm(:,7:8) = real(biome_mmte(:,7:8).^(1/1.3333)); %TELTL (1+1+2)/3 or (1.5+1+1.5)/3

q(:,1) = biome_te(:,1);
q(:,2) = biome_te(:,4);
q(:,3) = biome_te(:,3);
q(:,4) = biome_te(:,2);

Q = array2table(q,'VariableNames',{'TEM','TELTL','TEHTL','TEATL'},...
    'RowNames',{'LC','ECCS','ECSS'});

B = array2table(biome_teeff,'VariableNames',{'TEeffM','TEeffATL','TEeffHTL',...
    'TEeffLTL'},'RowNames',{'LC','ECCS','ECSS'});

% MAKE A TABLE OF MIN, MEAN, MAX
% ROWS = TLS
% COLS = MIN, MEAN, MAX BY BIOME
M = array2table(biome_mm,'VariableNames',{'minM','maxM','minATL','maxATL',...
    'minHTL','maxHTL','minLTL','maxLTL'},'RowNames',{'LC','ECCS','ECSS'});

% MAX TOO HIGH BECAUSE CAN EAT ~100% LOSS TO HP, NEED PROD INSTEAD

%% save
writetable(Q,[fpath 'TE_biomesMean_Historic9095_All_fish03_' cfile '.csv']...
    ,'Delimiter',',','WriteRowNames',true);
writetable(B,[fpath 'TEeff_biomes_Historic9095_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(M,[fpath 'TE_biomesMinMax_Historic9095_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);

