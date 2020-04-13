% Correlation analysis of pdiff in PRODUCTION
% of diff fn groups / trophic levels

%Kwia et al
%The coefficient and r2 of global phytoplankton?zooplankton linear 
%regressions. All regressions are based on global annual biomass anomalies and 
%are significant at the p < 0.05 level

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
dpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';

load([fpath 'Prod_pdiffs_hist_fore_',harv,'_' cfile '.mat']);

%% correlation analysis

Pdiff(:,1) = pdiffN(ID);
Pdiff(:,2) = pdiffZ(ID);
Pdiff(:,3) = pdiffDet(ID);
Pdiff(:,4) = pdiffF(ID);
Pdiff(:,5) = pdiffP(ID);
Pdiff(:,6) = pdiffD(ID);
Pdiff(:,7) = pdiffB(ID);

%%
[R,P] = corrcoef(Pdiff);

types = {'NPP','Z','Det','F','P','D','B'};
Rtab = array2table(R,'VariableNames',types,'RowNames',types);
Ptab = array2table(P,'VariableNames',types,'RowNames',types);

writetable(Rtab,[fpath 'Corr_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Rtab,[dpath 'Corr_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);

writetable(Ptab,[fpath 'CorrP_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab,[dpath 'CorrP_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);

%%
figure
pcolor(100*Pdiff)
shading flat
colorbar
cmocean('balance')
caxis([-100 100])
print('-dpng',[pp 'Hist_Fore_',harv,'_global_pdiff_prod_pcolor.png'])


