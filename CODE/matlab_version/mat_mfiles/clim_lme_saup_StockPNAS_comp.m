% Compare POEM vs. SAUP
% to Stock et al PNAS model vs. SAUP

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

%Reconciling Fisheries Catch and Ocean Productivity
%***TEMPLATE FOR FEEDBACK, PENDING FINAL CHECKS***
%model: 4
%alpha: 0.14
%TE_p: 0.13
%TE_b: 0.40
%f_T: 0.74
%T_100,warm: 19.99
%All fluxes in g C m-2 yr-1, Temperature in degrees celsius
%cols = LME  NPP   MESOZP  FDET   TLeq     T  modcatch SAUcatch
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'])

% POEM
%tab  = MT/km2
%tab2 = gC/m2
%cols = 'LME','ZP','Det','Bent','T','Pcatch','Scatch'
load([dpath 'LME_prod_catch_SAUP_' cfile '.mat'],'tab','tab2')

% Temps w/color
load([cpath 'LME_clim_temp_zoop_det.mat']);
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

%% Temp
figure(1)
subplot(2,3,1)
plot(-2:30,-2:30,'--k'); hold on;
plot(StockPNAS(:,6),tab(:,5),'.k','MarkerSize',20)
axis([-2 30 -2 30])
xlabel('Stock PNAS')
ylabel('POEM')
title('LME temp')

% Zoop
subplot(2,3,2)
plot(0:40,0:40,'--k'); hold on;
plot(StockPNAS(:,3),tab(:,2),'.k','MarkerSize',20)
xlabel('Stock PNAS')
ylabel('POEM')
title('LME mesozoop')

% Det
subplot(2,3,3)
plot(0:50,0:50,'--k'); hold on;
plot(StockPNAS(:,4),tab(:,3),'.k','MarkerSize',20)
axis([0 50 0 50])
xlabel('Stock PNAS')
ylabel('POEM')
title('LME bottom detritus')

% Modeled catch
subplot(2,3,4)
plot(0:0.1:1.2,0:0.1:1.2,'--k'); hold on;
plot(StockPNAS(:,7),tab2(:,6),'.k','MarkerSize',20)
axis([0 1.2 0 1.2])
xlabel('Stock PNAS')
ylabel('POEM')
title('LME model catch')

% SAUP catch
subplot(2,3,5)
plot(0:0.1:1.2,0:0.1:1.2,'--k'); hold on;
plot(StockPNAS(:,8),tab2(:,7),'.k','MarkerSize',20)
axis([0 1.2 0 1.2])
xlabel('Stock PNAS')
ylabel('POEM')
title('LME SAU catch')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_Stock_PNAS_all_comps.png'])

%% Comparison stats
%r
rall=corr(StockPNAS(:,7),tab2(:,6));
rall2=corr(StockPNAS(:,8),tab2(:,7));

%root mean square error
o=StockPNAS(:,7);
p=tab2(:,6);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=StockPNAS(:,8);
p=tab2(:,7);
n = length(o);
num=nansum((p-o).^2);
rmse2 = sqrt(num/n);

%Fmed
Fall=10^(median(StockPNAS(:,7)-tab2(:,6)));
Fall2=10^(median(StockPNAS(:,8)-tab2(:,7)));

% Table
fish_stat(1,1) = rall;
fish_stat(2,1) = rmse;
fish_stat(3,1) = Fall;
fish_stat(1,2) = rall2;
fish_stat(2,2) = rmse2;
fish_stat(3,2) = Fall2;

Fstat = array2table(fish_stat,'RowNames',{'r','RMSE','Fmed'},...
    'VariableNames',{'Models','SAU'});
writetable(Fstat,[dpath 'LME_Stock_PNAS_stats_' cfile '.csv'],...
    'Delimiter',',','WriteRowNames',true)
save([dpath 'LME_Stock_PNAS_stats_' cfile '.mat'],'fish_stat')

%% Figures
x=0:0.1:1.1;
% Correlation
figure(10)
plot(x,x,'--k');hold on;
for i=1:length(notLELC)
    lme=notLELC(i);
    plot(StockPNAS(i,7),tab2(i,6),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.05,0.95,['r = ' num2str(rall)])
text(0.05,0.9,['RMSE = ' num2str(rmse)])
text(0.05,0.85,['Fmed = ' num2str(Fall)])
axis([0 1.1 0 1.1])
xlabel('Stock model predict')
ylabel('POEM Climatology')
title('LME model catch')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_Stock_PNAS_catch_comp.png'])

% SAUP catch
figure(11)
plot(x,x,'--k');hold on;
for i=1:length(notLELC)
    lme=notLELC(i);
    plot(StockPNAS(i,8),tab2(i,7),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.05,0.95,['r = ' num2str(rall2)])
text(0.05,0.9,['RMSE = ' num2str(rmse2)])
text(0.05,0.85,['Fmed = ' num2str(Fall2)])
axis([0 1.1 0 1.1])
xlabel('Stock PNAS')
ylabel('POEM analyses')
title('LME SAU catch')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_Stock_PNAS_SAU_comp.png'])

