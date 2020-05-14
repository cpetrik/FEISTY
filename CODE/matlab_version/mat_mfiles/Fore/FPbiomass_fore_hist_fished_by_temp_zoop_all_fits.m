% Visualize biomass of F and P as fn of zoop
% at 2 different temps (or bins)
% ESM2M Hindcast of 1951-2000 and Forecast of 2051-2100

clear all
close all


%% FEISTY Output
% Original parameters
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];
harv = 'All_fish03';

fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
fpath2=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];

%load([fpath 'FP_zoop_temp_hist_fore_fits_' cfile '.mat'],'opred');
load([fpath2 'FP_zoop_temp_hist_fore_fits_' cfile '.mat'],'opred');

% Parameter ensemble, diff ks
efile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([efile 'FP_zoop_temp_hist_fore_fits_ens.mat'],'epred');

% Parameter enseble, same ks
nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'FP_zoop_temp_hist_fore_fits_ens_samek.mat'],'spred');

%% 10 = 7.5-12.5C, 20 = 22.5-27.5C
% opred(:,1) = pF10;
% opred(:,2) = pP10;
% opred(:,3) = pF20;
% opred(:,4) = pP20;

zoo = [0:25]';
X=[zoo; flipud(zoo)];

epred2 = cat(3,epred,opred);

% Cone of uncert 
emF10 = squeeze(mean(epred2(:,1,:),3));
emP10 = squeeze(mean(epred2(:,2,:),3));
emF20 = squeeze(mean(epred2(:,3,:),3));
emP20 = squeeze(mean(epred2(:,4,:),3));

esF10 = squeeze(std(epred2(:,1,:),0,3));
esP10 = squeeze(std(epred2(:,2,:),0,3));
esF20 = squeeze(std(epred2(:,3,:),0,3));
esP20 = squeeze(std(epred2(:,4,:),0,3));

%create y values for out and then back
%+/- 1 stdev
Ef10=[emF10+esF10; flipud(emF10-esF10)]; 
Ep10=[emP10+esP10; flipud(emP10-esP10)];
Ef20=[emF20+esF20; flipud(emF20-esF20)]; 
Ep20=[emP20+esP20; flipud(emP20-esP20)];


% Cone of uncert 
smF10 = squeeze(mean(spred(:,1,:),3));
smP10 = squeeze(mean(spred(:,2,:),3));
smF20 = squeeze(mean(spred(:,3,:),3));
smP20 = squeeze(mean(spred(:,4,:),3));

ssF10 = squeeze(std(spred(:,1,:),0,3));
ssP10 = squeeze(std(spred(:,2,:),0,3));
ssF20 = squeeze(std(spred(:,3,:),0,3));
ssP20 = squeeze(std(spred(:,4,:),0,3));

%create y values for out and then back
%+/- 1 stdev
Sf10=[smF10+ssF10; flipud(smF10-ssF10)]; 
Sp10=[smP10+ssP10; flipud(smP10-ssP10)];
Sf20=[smF20+ssF20; flipud(smF20-ssF20)]; 
Sp20=[smP20+ssP20; flipud(smP20-ssP20)];

%% same colors
figure(1)
subplot(2,2,1)
plot(zoo,opred(:,2),'b','LineWidth',2); hold on;
plot(zoo,opred(:,1),'r','LineWidth',2); hold on;
plot(zoo,squeeze(epred(:,2,:)),'b'); hold on;
plot(zoo,squeeze(epred(:,1,:)),'r'); hold on;
title('10^oC')
ylabel('Unequal temp-dep')
axis([0 25 0 20])
legend('P','F')
legend('location','northwest')

subplot(2,2,2)
plot(zoo,opred(:,4),'b','LineWidth',2); hold on;
plot(zoo,opred(:,3),'r','LineWidth',2); hold on;
plot(zoo,squeeze(epred(:,4,:)),'b'); hold on;
plot(zoo,squeeze(epred(:,3,:)),'r'); hold on;
title('25^oC')
axis([0 25 0 20])

subplot(2,2,3)
plot(zoo,squeeze(spred(:,2,:)),'b'); hold on;
plot(zoo,squeeze(spred(:,1,:)),'r'); hold on;
axis([0 25 0 20])
ylabel('Equal temp-dep')

subplot(2,2,4)
plot(zoo,squeeze(spred(:,4,:)),'b'); hold on;
plot(zoo,squeeze(spred(:,3,:)),'r'); hold on;
text(-4,-4,'Zooplankton biomass (g m^-^2)','HorizontalAlignment','center','FontWeight','bold')
text(-39,23,'Fish biomass (g m^-^2)','HorizontalAlignment','center','Rotation',90,'FontWeight','bold')
axis([0 25 0 20])
print('-dpng',[pp 'Hist_Fore_',harv,'_FP_zoop_temp_fits_lines_v1.png'])

%% lighter colors for ensemble
figure(2)
subplot(2,2,1)
plot(zoo,opred(:,2),'b','LineWidth',2); hold on;
plot(zoo,opred(:,1),'r','LineWidth',2); hold on;
plot(zoo,squeeze(epred(:,2,:)),'color',[0/255 240/255 240/255]); hold on;
plot(zoo,squeeze(epred(:,1,:)),'color',[255/255 192/255 203/255]); hold on;
plot(zoo,opred(:,2),'b','LineWidth',2); hold on;
plot(zoo,opred(:,1),'r','LineWidth',2); hold on;
title('10^oC')
ylabel('Unequal temp-dep')
axis([0 25 0 20])
legend('P','F')
legend('location','northwest')

subplot(2,2,2)
plot(zoo,squeeze(epred(:,4,:)),'color',[0/255 240/255 240/255]); hold on;
plot(zoo,squeeze(epred(:,3,:)),'color',[255/255 192/255 203/255]); hold on;
plot(zoo,opred(:,4),'b','LineWidth',2); hold on;
plot(zoo,opred(:,3),'r','LineWidth',2); hold on;
title('25^oC')
axis([0 25 0 20])

subplot(2,2,3)
plot(zoo,squeeze(spred(:,2,:)),'color',[0/255 240/255 240/255]); hold on;
plot(zoo,squeeze(spred(:,1,:)),'color',[255/255 192/255 203/255]); hold on;
axis([0 25 0 20])
ylabel('Equal temp-dep')

subplot(2,2,4)
plot(zoo,squeeze(spred(:,4,:)),'color',[0/255 240/255 240/255]); hold on;
plot(zoo,squeeze(spred(:,3,:)),'color',[255/255 192/255 203/255]); hold on;
text(-4,-4,'Zooplankton biomass (g m^-^2)','HorizontalAlignment','center','FontWeight','bold')
text(-39,23,'Fish biomass (g m^-^2)','HorizontalAlignment','center','Rotation',90,'FontWeight','bold')
axis([0 25 0 20])
print('-dpng',[pp 'Hist_Fore_',harv,'_FP_zoop_temp_fits_lines_v2.png'])

%% +/- std dev
figure(3)
subplot(2,2,1)
plot(zoo,emP10,'b','LineWidth',2); hold on;
plot(zoo,emF10,'r','LineWidth',2); hold on;
fill(X,Ep10,'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,Ef10,'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
title('10^oC')
ylabel('Unequal temp-dep')
axis([0 25 0 20])
legend('P','F')
legend('location','northwest')
%
subplot(2,2,2)
plot(zoo,emP20,'b','LineWidth',2); hold on;
plot(zoo,emF20,'r','LineWidth',2); hold on;
fill(X,Ep20,'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,Ef20,'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
title('25^oC')
axis([0 25 0 20])

subplot(2,2,3)
plot(zoo,smP10,'b','LineWidth',2); hold on;
plot(zoo,smF10,'r','LineWidth',2); hold on;
fill(X,Sp10,'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,Sf10,'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
axis([0 25 0 20])
ylabel('Equal temp-dep')

subplot(2,2,4)
plot(zoo,smP20,'b','LineWidth',2); hold on;
plot(zoo,smF20,'r','LineWidth',2); hold on;
fill(X,Sp20,'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,Sf20,'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
text(-4,-4,'Zooplankton biomass (g m^-^2)','HorizontalAlignment','center','FontWeight','bold')
text(-39,23,'Fish biomass (g m^-^2)','HorizontalAlignment','center','Rotation',90,'FontWeight','bold')
axis([0 25 0 20])
print('-dpng',[pp 'Hist_Fore_',harv,'_FP_zoop_temp_fits_cone.png'])




