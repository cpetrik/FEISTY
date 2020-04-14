% Correlation analysis of pdiff in PRODUCTION
% of diff fn groups / trophic levels

%Kwia et al
%The coefficient and r2 of global phytoplankton-zooplankton linear 
%regressions. All regressions are based on global annual biomass anomalies and 
%are significant at the p < 0.05 level

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); 
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

%max(pdiffF(ID)) = 4.0130e+07
%max(pdiffP(ID)) = 2.0409e+05
Pdiff(Pdiff > 2) = nan;         %F & P corrs sensitive to this choice!

%% Correlations
[R,P] = corrcoef(Pdiff,'Rows','pairwise');
[Rp,Pp] = corr(Pdiff,'Rows','pairwise','Type','Pearson');
[Rk,Pk] = corr(Pdiff,'Rows','pairwise','Type','Kendall');
[Rs,Ps] = corr(Pdiff,'Rows','pairwise','Type','Spearman');

%%
types = {'NPP','Z','Det','F','P','D','B'};
Rtab = array2table(R,'VariableNames',types,'RowNames',types);
Ptab = array2table(P,'VariableNames',types,'RowNames',types);

Rktab = array2table(Rk,'VariableNames',types,'RowNames',types);
Pktab = array2table(Pk,'VariableNames',types,'RowNames',types);

Rstab = array2table(Rs,'VariableNames',types,'RowNames',types);
Pstab = array2table(Ps,'VariableNames',types,'RowNames',types);

writetable(Rtab,[fpath 'Corr_Pearson_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Rtab,[dpath 'Corr_Pearson_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab,[fpath 'CorrP_Pearson_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab,[dpath 'CorrP_Pearson_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);

writetable(Rktab,[fpath 'Corr_Kendall_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Rktab,[dpath 'Corr_Kendall_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Pktab,[fpath 'CorrP_Kendall_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Pktab,[dpath 'CorrP_Kendall_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);

writetable(Rstab,[fpath 'Corr_Spearman_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Rstab,[dpath 'Corr_Spearman_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Pstab,[fpath 'CorrP_Spearman_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Pstab,[dpath 'CorrP_Spearman_Fore_Hist_All_fish03_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true);

%% Linear regressions
%NPP
nz = fitlm(Pdiff(:,1),Pdiff(:,2))      %Z = -0.02292 + 1.7699N
ne = fitlm(Pdiff(:,1),Pdiff(:,3))      %E = -0.27895 + 0.92944N
nf = fitlm(Pdiff(:,1),Pdiff(:,4))      %F = -0.00918 + 0.60175N
np = fitlm(Pdiff(:,1),Pdiff(:,5))      %P = -0.56157 + 0.38446N
nd = fitlm(Pdiff(:,1),Pdiff(:,6))      %D = -0.38528 + 1.1334N
nb = fitlm(Pdiff(:,1),Pdiff(:,7))      %B = -0.34664 + 1.2279N

%Zoop
ze = fitlm(Pdiff(:,2),Pdiff(:,3))   %E = -0.27175  + 0.4343Z
zf = fitlm(Pdiff(:,2),Pdiff(:,4))   %F =  0.023824 + 0.77374Z
zp = fitlm(Pdiff(:,2),Pdiff(:,5))   %P = -0.5453   + 0.41157Z
zd = fitlm(Pdiff(:,2),Pdiff(:,6))   %D = -0.37559  + 0.54559Z
zb = fitlm(Pdiff(:,2),Pdiff(:,7))   %B = -0.33623  + 0.57741Z

%Det
ef = fitlm(Pdiff(:,3),Pdiff(:,4))   %F =  0.078142 + 0.32933E
ep = fitlm(Pdiff(:,3),Pdiff(:,5))   %P = -0.45573  + 0.38203E
ed = fitlm(Pdiff(:,3),Pdiff(:,6))   %D = -0.043686 + 1.2243E
eb = fitlm(Pdiff(:,3),Pdiff(:,7))   %B =  0.019228 + 1.3119E

%F
fp = fitlm(Pdiff(:,4),Pdiff(:,5))   %P = -0.55968 + 0.19623F
fd = fitlm(Pdiff(:,4),Pdiff(:,6))   %D = -0.40249 + 0.036706F
fb = fitlm(Pdiff(:,4),Pdiff(:,7))   %B = -0.36548 + 0.020894F

%P
pd = fitlm(Pdiff(:,5),Pdiff(:,6))   %D = -0.38935 + 0.024974P
pb = fitlm(Pdiff(:,5),Pdiff(:,7))   %B = -0.35515 + 0.020091P

%D
db = fitlm(Pdiff(:,6),Pdiff(:,7))   %B = 0.047866 + 1.0265D

%B
bd = fitlm(Pdiff(:,7),Pdiff(:,6))   %D = -0.13381 + 0.73639B

%%
Pdiff(:,8) = nan(size(Pdiff(:,7)));
figure
pcolor(100*Pdiff)
shading flat
colorbar
cmocean('balance')
caxis([-100 100])
set(gca,'XTick',[1.5:7.5],'XTickLabel',types)
print('-dpng',[pp 'Hist_Fore_',harv,'_global_pdiff_prod_pcolor.png'])


