% Compare J&C15 to POEM and COBALT
% COBALT data from Stock et al PNAS 

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

%cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
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
CPP = StockPNAS(:,2)./365;

% POEM
load([dpath 'LME_clim_fished_' harv '_' cfile '.mat'])
load([dpath 'LME_prod_catch_SAUP_' cfile '.mat'],'tab','tab2')
%lme_mbio = mean biomass in g WW
%lme_sbio = total biomass in g WW
lme_mbio_MT = lme_mbio.*1e-6;
lme_sbio_MT = lme_sbio.*1e-6;
%sum all fishes
lme_tmbio_MT = sum(lme_mbio_MT(:,1:8),2);
lme_tsbio_MT = sum(lme_sbio_MT(:,1:8),2);
lme_tmbio_logMT = log10(lme_tmbio_MT);
lme_tsbio_logMT = log10(lme_tsbio_MT);
[sort_mbio,plme_rank_mbio] = sort(lme_tmbio_logMT,'descend');

% J&C15 data
%log10 total consumer biomass in 10^6 tonnes (ranges -1 to 3)
%10^6 tonnes = 1 MT = 10^6 g
%PP in g C m-2 d-1
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'J&C15_all.mat'])

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
subplot(1,2,1)
plot(-2:30,-2:30,'--k'); hold on;
plot(StockPNAS(:,6),JCT(notLELC),'.k','MarkerSize',20)
axis([-2 30 -2 30])
axis equal
xlabel('COBALT (Stock PNAS)')
ylabel('J&C15')
title('LME temp')

subplot(1,2,2)
plot(-2:30,-2:30,'--k'); hold on;
plot(tab(:,5),JCT(notLELC),'.k','MarkerSize',20)
axis([-2 30 -2 30])
axis equal
xlabel('COBALT (POEM)')
ylabel('J&C15')
title('LME temp')
print('-dpng',[ppath 'Clim_',harv,'_JC15_temp_comp.png'])

%% NPP
JCPP2 = JCPP(notLELC);
figure(2)
plot(0:0.2:1.8,0:0.2:1.8,'--k'); hold on;
for i=1:length(notLELC)
    lme=notLELC(i);
    plot(CPP(i),JCPP2(i),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
xlabel('COBALT (Stock PNAS)')
ylabel('J&C15')
title('LME NPP (gC m^-^2 d^-^1)')
print('-dpng','/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/JC15_NPP_comp.png')

figure(4)
plot(0:0.2:1.8,0:0.2:1.8,'--k'); hold on;
for i=1:length(notLELC)
    lme=notLELC(i);
    plot(CPP(i),JCPP2(i),'ok','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(CPP(i),JCPP2(i),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
xlabel('COBALT (Stock PNAS)')
ylabel('J&C15')
title('LME NPP (gC m^-^2 d^-^1)')
print('-dpng','/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/JC15_NPP_compID.png')

%% Modeled biomass rank
figure(3)
plot(1:66,1:66,'--k'); hold on;
for i=1:66
    plot(plme_rank_mbio(i),JCrank(i),'.k','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
end
ylabel('J&C15')
xlabel('POEM')
title('LME ranking high to low biomass')
print('-dpng',[ppath 'Clim_',harv,'_JC15_rank_comp.png'])

%% Stats in R
%r
rall=corr(StockPNAS(:,6),JCT(notLELC)); %0.98
rall2=corr(tab(:,5),JCT(notLELC));      %0.98
rall3=corr(CPP,JCPP(notLELC));          %0.54

%root mean square error
o=StockPNAS(:,6);
p=JCT(notLELC);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);     %2.7319

o=tab(:,5);
p=JCT(notLELC);
n = length(o);
num=nansum((p-o).^2);
rmse2 = sqrt(num/n);    %2.7702

o=CPP;
p=JCPP(notLELC);
n = length(o);
num=nansum((p-o).^2);
rmse3 = sqrt(num/n);    %0.3597

%Fmed
Fall=10^(median(StockPNAS(:,6)-JCT(notLELC)));  %0.0135
Fall2=10^(median(tab(:,5)-JCT(notLELC)));       %0.0125
Fall3=10^(median(CPP-JCPP(notLELC)));           %0.8486

% Table
fish_stat(1,1) = rall;
fish_stat(2,1) = rmse;
fish_stat(3,1) = Fall;
fish_stat(1,2) = rall2;
fish_stat(2,2) = rmse2;
fish_stat(3,2) = Fall2;
fish_stat(1,3) = rall3;
fish_stat(2,3) = rmse3;
fish_stat(3,3) = Fall3;

Fstat = array2table(fish_stat,'RowNames',{'r','RMSE','Fmed'},...
    'VariableNames',{'StockT','POEMT','StockNPP'});
writetable(Fstat,[cpath 'LME_Stock_PNAS_JC15_stats.csv'],...
    'Delimiter',',','WriteRowNames',true)
save([cpath 'LME_Stock_PNAS_JC15_stats.mat'],'fish_stat')

Rstat = table(LME,JCT,JCPP,JCrank,lme_ptemp(:,1),plme_rank_mbio,...
    'VariableNames',{'LME','JCT','JCPP','JCrank','PT','Prank'});
Lstat = table(notLELC,JCT(notLELC),JCPP2,JCrank(notLELC),CPP,StockPNAS(:,6),...
    'VariableNames',{'LME','JCT','JCPP','JCrank','CNPP','CT'});
writetable(Rstat,[dpath 'LME_JC15_values_',cfile,'.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Lstat,[cpath 'LME_Stock_PNAS_JC15_values.csv'],...
    'Delimiter',',','WriteRowNames',true)



