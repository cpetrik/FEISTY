% 1961-2010 obsclim
% Observed effort
% Compare diff versions

clear 
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];


pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP/Phase3a/';
ppath = [pp cfile '/OneDeg/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%% v1
mod = 'obsclim_All_fishobs_'; %'obsclim_All_fishobs_v3_';
load([fpath 'LME_Hist_',mod,cfile,'.mat']);

lme_area_km2 = lme_area * 1e-6;

%mcatch = mf, mp, md, lp, ld
%MT/km2
Alme_fcatch_all = nansum(lme_mcatch,3) * 1e-6 ./ lme_area_km2;
Flme_fcatch_all = (lme_mcatch(:,:,1)) * 1e-6 ./ lme_area_km2;
Plme_fcatch_all = (lme_mcatch(:,:,2)+lme_mcatch(:,:,4)) * 1e-6 ./ lme_area_km2;
Dlme_fcatch_all = (lme_mcatch(:,:,3)+lme_mcatch(:,:,5)) * 1e-6 ./ lme_area_km2;

Fv1 = sum(Flme_fcatch_all);
Pv1 = sum(Plme_fcatch_all);
Dv1 = sum(Dlme_fcatch_all);
Av1 = sum(Alme_fcatch_all);

%% v1.2
mod = 'obsclim_All_fishobs_v1.2_';
load([fpath 'LME_Hist_',mod,cfile,'.mat']);

lme_area_km2 = lme_area * 1e-6;

%mcatch = mf, mp, md, lp, ld
%MT/km2
Alme_fcatch_all = nansum(lme_mcatch,3) * 1e-6 ./ lme_area_km2;
Flme_fcatch_all = (lme_mcatch(:,:,1)) * 1e-6 ./ lme_area_km2;
Plme_fcatch_all = (lme_mcatch(:,:,2)+lme_mcatch(:,:,4)) * 1e-6 ./ lme_area_km2;
Dlme_fcatch_all = (lme_mcatch(:,:,3)+lme_mcatch(:,:,5)) * 1e-6 ./ lme_area_km2;

Fv12 = sum(Flme_fcatch_all);
Pv12 = sum(Plme_fcatch_all);
Dv12 = sum(Dlme_fcatch_all);
Av12 = sum(Alme_fcatch_all);

%% v2
mod = 'obsclim_All_fishobs_v2_';
load([fpath 'LME_Hist_',mod,cfile,'.mat']);

lme_area_km2 = lme_area * 1e-6;

%mcatch = mf, mp, md, lp, ld
%MT/km2
Alme_fcatch_all = nansum(lme_mcatch,3) * 1e-6 ./ lme_area_km2;
Flme_fcatch_all = (lme_mcatch(:,:,1)) * 1e-6 ./ lme_area_km2;
Plme_fcatch_all = (lme_mcatch(:,:,2)+lme_mcatch(:,:,4)) * 1e-6 ./ lme_area_km2;
Dlme_fcatch_all = (lme_mcatch(:,:,3)+lme_mcatch(:,:,5)) * 1e-6 ./ lme_area_km2;

Fv2 = sum(Flme_fcatch_all);
Pv2 = sum(Plme_fcatch_all);
Dv2 = sum(Dlme_fcatch_all);
Av2 = sum(Alme_fcatch_all);

%% v3
mod = 'obsclim_All_fishobs_v3_';
load([fpath 'LME_Hist_',mod,cfile,'.mat']);

lme_area_km2 = lme_area * 1e-6;

%mcatch = mf, mp, md, lp, ld
%MT/km2
Alme_fcatch_all = nansum(lme_mcatch,3) * 1e-6 ./ lme_area_km2;
Flme_fcatch_all = (lme_mcatch(:,:,1)) * 1e-6 ./ lme_area_km2;
Plme_fcatch_all = (lme_mcatch(:,:,2)+lme_mcatch(:,:,4)) * 1e-6 ./ lme_area_km2;
Dlme_fcatch_all = (lme_mcatch(:,:,3)+lme_mcatch(:,:,5)) * 1e-6 ./ lme_area_km2;

Fv3 = sum(Flme_fcatch_all);
Pv3 = sum(Plme_fcatch_all);
Dv3 = sum(Dlme_fcatch_all);
Av3 = sum(Alme_fcatch_all);

%% SAUP catch
spath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/SAUP/';
%use weighted catches
load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
    'Flme_wcatch','Dlme_wcatch','Plme_wcatch');

Flme_scatch_all = nansum(Flme_wcatch,3);
Plme_scatch_all = nansum(Plme_wcatch,3);
Dlme_scatch_all = nansum(Dlme_wcatch,3);

%% SAUP 1950-2010
% Phase3a 1961-2010
years = 1961:2010;
%years = 1961:2005;
id = find(yr>=years(1) & yr<=years(end));

slme_scatch_all = lme_catch(id,:);
Flme_scatch_all = Flme_scatch_all(id,:);
Plme_scatch_all = Plme_scatch_all(id,:);
Dlme_scatch_all = Dlme_scatch_all(id,:);
slme_scatch_all = slme_scatch_all' ./ lme_area_km2;
Flme_scatch_all = Flme_scatch_all' ./ lme_area_km2;
Plme_scatch_all = Plme_scatch_all' ./ lme_area_km2;
Dlme_scatch_all = Dlme_scatch_all' ./ lme_area_km2;

SF = sum(Flme_scatch_all);
SP = sum(Plme_scatch_all);
SD = sum(Dlme_scatch_all);
SA = sum(slme_scatch_all);

%% Time series Plots 
% Global
figure(7)
subplot(2,2,1)
plot(years,SF,'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,Fv1,'k','LineWidth',2)
plot(years,Fv12,'r','LineWidth',2)
plot(years,Fv2,'b','LineWidth',2)
plot(years,Fv3,'color',[0 0.6 0],'LineWidth',2)
legend('SAU','v1','v1.2','v2','v3')
legend('location','northwest')
ylabel('FEISTY catch (MT)')
title('Forage')
%ylim([0 60])

subplot(2,2,2)
plot(years,SP,'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,Pv1,'k','LineWidth',2)
plot(years,Pv12,'r','LineWidth',2)
plot(years,Pv2,'b','LineWidth',2)
plot(years,Pv3,'color',[0 0.6 0],'LineWidth',2)
title('Large pelagic')
%ylim([10 60])

subplot(2,2,3)
plot(years,SD,'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,Dv1,'k','LineWidth',2)
plot(years,Dv12,'r','LineWidth',2)
plot(years,Dv2,'b','LineWidth',2)
plot(years,Dv3,'color',[0 0.6 0],'LineWidth',2)
xlabel('year')
%ylabel('FEISTY catch (MT km^-^2)')
ylabel('FEISTY catch (MT)')
title('Demersal')
%ylim([25 70])

subplot(2,2,4)
plot(years,SA,'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,Av1,'k','LineWidth',2)
plot(years,Av12,'r','LineWidth',2)
plot(years,Av2,'b','LineWidth',2)
plot(years,Av3,'color',[0 0.6 0],'LineWidth',2)
xlabel('year')
title('All')
%ylim([40 180])
stamp(mod)
print('-dpng',[ppath 'Hist_obsclim_All_fishobs_comp_onedeg_sumLME_ts.png'])

%% log10 Global
figure(8)
subplot(2,2,1)
plot(years,log10(SF),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,log10(Fv1),'k','LineWidth',2)
plot(years,log10(Fv12),'r','LineWidth',2)
plot(years,log10(Fv2),'b','LineWidth',2)
plot(years,log10(Fv3),'color',[0 0.6 0],'LineWidth',2)
legend('SAU','v1','v1.2','v2','v3')
legend('location','northwest')
ylabel('log_1_0 FEISTY catch (MT km)')
title('Forage')
%ylim([0.5 2])

subplot(2,2,2)
plot(years,log10(SP),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,log10(Pv1),'k','LineWidth',2)
plot(years,log10(Pv12),'r','LineWidth',2)
plot(years,log10(Pv2),'b','LineWidth',2)
plot(years,log10(Pv3),'color',[0 0.6 0],'LineWidth',2)
title('Large pelagic')
%ylim([1 1.8])

subplot(2,2,3)
plot(years,log10(SD),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,log10(Dv1),'k','LineWidth',2)
plot(years,log10(Dv12),'r','LineWidth',2)
plot(years,log10(Dv2),'b','LineWidth',2)
plot(years,log10(Dv3),'color',[0 0.6 0],'LineWidth',2)
ylabel('year')
%ylabel('log_1_0 FEISTY catch (MT km^-^2)')
ylabel('log_1_0 FEISTY catch (MT km)')
title('Demersal')
%ylim([1.45 1.85])

subplot(2,2,4)
plot(years,log10(SA),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,log10(Av1),'k','LineWidth',2)
plot(years,log10(Av12),'r','LineWidth',2)
plot(years,log10(Av2),'b','LineWidth',2)
plot(years,log10(Av3),'color',[0 0.6 0],'LineWidth',2)
xlabel('year')
title('All')
%ylim([1.6 2.3])
stamp(mod)
print('-dpng',[ppath 'Hist_obsclim_All_fishobs_comp_onedeg_sumLME_ts_log10.png'])

