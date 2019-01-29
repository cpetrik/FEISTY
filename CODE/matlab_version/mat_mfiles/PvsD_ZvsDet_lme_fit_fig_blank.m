% Plot mean biomasses of zoop to compare to POEM results
% COMPARE BY LME, COLOR-CODE BY TEMP, Plot GAM fit
% Uses 3 Knots and Probit link fn

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

gpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
cpath = [Pdrpbx 'Princeton/POEM_other/cobalt_data/'];
pp = [Pdrpbx 'Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([gpath 'esm26_area_1deg.mat']);
load([gpath 'LME_clim_temp_zoop_det_npp.mat']);
load([gpath 'LME_depth_area.mat'],'lme_depth','lme_shal_frac');

%% POEM results
frate = 0.3;
tfish = num2str(100+int64(10*frate));
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = ['All_fish',tfish(2:end)];
tharv = ['Harvest all fish ' num2str(frate)];
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'ZlDet_PD_gam_fit_extrap.mat']);

%% Scatter plot
% Plot GAM fit
figure(1)
plot(ZlDet,PDfit,'k'); hold on;
plot(ZlDet,PDfit+2*PDse,'--k'); hold on;
plot(ZlDet,PDfit-2*PDse,'--k'); hold on;
%title('Color = LME mean temp')
xlabel('F_{photic}:F_{benthic}')
ylabel('Fraction of pelagic fish')
axis([-1 0.75 0 1])
print('-dpng',[ppath 'lme_scatter_ZlDet_PD_GAMfit_BW.png'])

%% Scatter plot
PDse2 = PDse;
id = (PDse > mean(PDse));
PDse2(id) = mean(PDse);
% Mean SE
figure(2)
plot(ZlDet,PDfit,'k','LineWidth',2); hold on;
plot(ZlDet,PDfit+2*PDse2,'--k'); hold on;
plot(ZlDet,PDfit-2*PDse2,'--k'); hold on;
%title('Color = LME mean temp')
xlabel('F_{pelagic}:F_{benthic}')
ylabel('Fraction of pelagic fish')
%axis([-1 0.75 0 1])
set(gca,'XTick',[-1:1],'XTickLabel',[1,10,100],'FontSize',14)
print('-dpng',[ppath 'lme_scatter_ZlDet_PD_GAMfit_BW_meanSE.png'])

%% Fake curve
N=1:0.1:3;
Frac = 0.1*N.^3 ./ (1+0.1*N.^3);

figure
plot(N,Frac,'k'); hold on;


