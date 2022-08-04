% test diff

clear all
close all

%%
load('/Volumes/MIP/Fish-MIP/Phase3/fishing/grid_mortality_guilds/CORE_mortality_all_ID_annual_tempSc.mat',...
    'fmD','fmF','fmP');
sfD = fmD;
sfF = fmF;
sfP = fmP;
clear fmD fmF fmP

% load(['/Users/cpetrik/Dropbox/Princeton/FEISTY_other/fishing_ms_ideas/fishing_effort_impl/grid_mortality_guilds/'...
%     'CORE_mortality_all_ID_annual_tempSc.mat'],...
%     'year','fmD','fmF','fmP');
load('/Volumes/MIP/Fish-MIP/Phase3/fishing/grid_mortality_guilds/CORE_mortality_all_ID_annual.mat',...
    'fmD','fmF','fmP','year');

%% Plots in time
y = year;
F1 = mean(sfF);
P1 = mean(sfP);
D1 = mean(sfD);
F2 = mean(fmF);
P2 = mean(fmP);
D2 = mean(fmD);

% All size classes of all
figure(1)
plot(y,(F1),'r','Linewidth',2); hold on;
plot(y,(P1),'b','Linewidth',2); hold on;
plot(y,(D1),'k','Linewidth',2); hold on;
plot(y,(F2),'--m','Linewidth',2); hold on;
plot(y,(P2),'--c','Linewidth',2); hold on;
plot(y,(D2),'--g','Linewidth',2); hold on;
% legend('F','P','D')
% legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([0.01 0.11])
xlabel('Year')
ylabel('mean fishing mortality (y^-^1)')
title(['CORE'])
%print('-dpng',[ppath 'CORE_ts_fishing_mort',harv,'_all_types.png'])