% Visualize time series output of FEISTY ensembles
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Saved as mat files
% Ensemble mid6, temp3
% Production & fisheries yield together on subplot

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];

efile = 'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050';
epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',efile,'/'];
epath2 = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/',efile,'/'];

%% Yield
load([epath 'Hist_Fore_All_fish03_ensem6_mid_temp3_ts_yield.mat'])
test=find(y>1950);
yid=test(1);

% types - diff
figure(1)
subplot(2,2,2)
plot(y,(rmF)*1e-6,'r','LineWidth',2); hold on;
plot(y,(rmP)*1e-6,'b','LineWidth',2); hold on;
plot(y,(rmD)*1e-6,'color',[0 0.6 0],'LineWidth',2); hold on;
fill(X,(Rf)*1e-6,'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Rp)*1e-6,'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Rd)*1e-6,'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
xlim([y(yid) y(end)])
ylim([-20 5])
text(y(yid),6.5,'B')
xlabel('Year')
ylabel('Yield (MT km^-^2 y^-^1) relative to 1951')
legend('F','P','D')
legend('location','southwest')

% adults - diff
figure(2)
subplot(2,2,2)
plot(y,(rmF)*1e-6,'r','LineWidth',2); hold on;
plot(y,(rmLP)*1e-6,'b','LineWidth',2); hold on;
plot(y,(rmLD)*1e-6,'color',[0 0.6 0],'LineWidth',2); hold on;
fill(X,(Rf)*1e-6,'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Rlp)*1e-6,'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Rld)*1e-6,'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
xlim([y(yid) y(end)])
ylim([-20 5])
text(y(yid),6.5,'B')
xlabel('Year')
ylabel('Yield (MT km^-^2 y^-^1) relative to 1951')
legend('F','P','D')
legend('location','southwest')

%% var in yield
hyr = (y>=1951 & y<2000);
fyr = (y>=2051 & y<2100);

vt(1,1) = mean(var(dtF(:,hyr),0,2)); %2.27         All * 1e+12
vt(1,2) = mean(var(dtF(:,fyr),0,2)); %1.47 F decr
vt(2,1) = mean(var(dtLP(:,hyr),0,2)); %1.77
vt(2,2) = mean(var(dtLP(:,fyr),0,2)); %2.67 P incr
vt(3,1) = mean(var(dtLD(:,hyr),0,2)); %0.26
vt(3,2) = mean(var(dtLD(:,fyr),0,2))  %0.83 D incr

%% Production
load([epath2 'Hist_Fore_All_fish03_ensem6_mid_temp3_ts_prod.mat'])
test=find(y>1950);
yid=test(1);

figure(1)
subplot(2,2,1)
fill(X,100*(Vf),'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,100*(Vp),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,100*(Vd),'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,100*(mpF),'r','LineWidth',2); hold on;
plot(y,100*(mpP),'b','LineWidth',2); hold on;
plot(y,100*(mpD),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
ylim([-35 15])
text(y(yid),18,'A')
xlabel('Year')
ylabel('% \Delta in production relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_yield_prod_types_all_ensem_mid6_temp3_cone_1std_yr_2plot.png'])

figure(2)
subplot(2,2,1)
fill(X,100*(Vf),'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,100*(Vp),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,100*(Vd),'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,100*(mpF),'r','LineWidth',2); hold on;
plot(y,100*(mpP),'b','LineWidth',2); hold on;
plot(y,100*(mpD),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
ylim([-35 15])
text(y(yid),18,'A')
xlabel('Year')
ylabel('% \Delta in production relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_yield_prod_types_adults_ensem_mid6_temp3_cone_1std_yr_2plot.png'])


%% var in prod
hyr = (y>=1951 & y<2000);
fyr = (y>=2051 & y<2100);

vt(4,1) = mean(var(dtF(:,hyr),0,2)); %0.21          All *1e-5
vt(4,2) = mean(var(dtF(:,fyr),0,2)); %0.21 F same
vt(5,1) = mean(var(dtP(:,hyr),0,2)); %0.02
vt(5,2) = mean(var(dtP(:,fyr),0,2)); %0.01 P decr
vt(6,1) = mean(var(dtD(:,hyr),0,2)); %0.0006
vt(6,2) = mean(var(dtD(:,fyr),0,2))  %0.0010 D incr

Stab = array2table(vt,'RowNames',{'Fyield','Pyield','Dyield','Fprod',...
    'Pprod','Dprod'},'VariableNames',{'Hist','Proj'});
writetable(Stab,[epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_tsVar_prod_yield.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Stab,[fpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_tsVar_prod_yield.csv'],...
    'Delimiter',',','WriteRowNames',true)


