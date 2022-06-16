% CORE-forced
% Observed effort

clear all
close all

%%
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/CORE/'];
ppath = [pp cfile '/CORE/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

harv = 'fished_obs';
tharv = 'Observed effort';

load([fpath 'LME_core_',harv,'_' cfile '.mat']);

years = 1961:2007;

%% FEISTY total catches
lme_area_km2 = lme_area * 1e-6;

%mcatch = mf, mp, md, lp, ld
%MT/km2
Alme_fcatch_all = nansum(lme_mcatch,3) * 1e-6 ./ lme_area_km2;
Flme_fcatch_all = (lme_mcatch(:,:,1)) * 1e-6 ./ lme_area_km2;
Plme_fcatch_all = (lme_mcatch(:,:,2)+lme_mcatch(:,:,4)) * 1e-6 ./ lme_area_km2;
Dlme_fcatch_all = (lme_mcatch(:,:,3)+lme_mcatch(:,:,5)) * 1e-6 ./ lme_area_km2;

%% SAUP catch
spath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/SAUP/';
%use weighted catches
load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
    'Flme_wcatch','Dlme_wcatch','Plme_wcatch');

Flme_scatch_all = nansum(Flme_wcatch,3);
Plme_scatch_all = nansum(Plme_wcatch,3);
Dlme_scatch_all = nansum(Dlme_wcatch,3);

%1961-2007 in SAUP 
id = find(yr>=years(1) & yr<=years(end));

slme_scatch_all = lme_catch(id,:);
Flme_scatch_all = Flme_scatch_all(id,:);
Plme_scatch_all = Plme_scatch_all(id,:);
Dlme_scatch_all = Dlme_scatch_all(id,:);
slme_scatch_all = slme_scatch_all' ./ lme_area_km2;
Flme_scatch_all = Flme_scatch_all' ./ lme_area_km2;
Plme_scatch_all = Plme_scatch_all' ./ lme_area_km2;
Dlme_scatch_all = Dlme_scatch_all' ./ lme_area_km2;

%% Scatter Plots 
% Global
figure(1)
subplot(2,2,1)
plot(Flme_scatch_all(:),Flme_fcatch_all(:),'.r','MarkerSize',20)
ylabel('FEISTY catch (MT km^-^2)')
title('Forage')

subplot(2,2,2)
plot(Plme_scatch_all(:),Plme_fcatch_all(:),'.b','MarkerSize',20)
title('Large pelagic')

subplot(2,2,3)
plot(Dlme_scatch_all(:),Dlme_fcatch_all(:),'.','color',[0 0.6 0],'MarkerSize',20)
xlabel('SAUP')
ylabel('FEISTY catch (MT km^-^2)')
title('Demersal')

subplot(2,2,4)
plot(slme_scatch_all(:),Alme_fcatch_all(:),'.k','MarkerSize',20)
xlabel('SAUP')
title('All')
print('-dpng',[ppath 'CORE_obs_fished_lme_scatter_all_yrs.png'])

%% log10 Global
figure(2)
subplot(2,2,1)
plot(log10(Flme_scatch_all(:)),log10(Flme_fcatch_all(:)),'.r','MarkerSize',20)
ylabel('log_1_0 FEISTY catch (MT km^-^2)')
title('Forage')

subplot(2,2,2)
plot(log10(Plme_scatch_all(:)),log10(Plme_fcatch_all(:)),'.b','MarkerSize',20)
title('Large pelagic')

subplot(2,2,3)
plot(log10(Dlme_scatch_all(:)),log10(Dlme_fcatch_all(:)),'.','color',[0 0.6 0],'MarkerSize',20)
xlabel('SAUP')
ylabel('log_1_0 FEISTY catch (MT km^-^2)')
title('Demersal')

subplot(2,2,4)
plot(log10(slme_scatch_all(:)),log10(Alme_fcatch_all(:)),'.k','MarkerSize',20)
xlabel('SAUP')
title('All')
print('-dpng',[ppath 'CORE_obs_fished_lme_scatter_all_yrs_log10.png'])

%% Scatter by LME
nrows=11;
ncols=6;
pos = subfigrid(nrows,ncols,[0.05 0.025 0.025 0.05],[0.72 0.75]);

%% all 
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(slme_scatch_all(i,:),Alme_fcatch_all(i,:),'.k','MarkerSize',10); hold on;
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'CORE_byLME_obs_fished_scatter_all_yrs_Allfish.png'])

%% F
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(Flme_scatch_all(i,:),Flme_fcatch_all(i,:),'.r','MarkerSize',10); hold on;
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'CORE_byLME_obs_fished_scatter_all_yrs_F.png'])

%% P
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(Plme_scatch_all(i,:),Plme_fcatch_all(i,:),'.b','MarkerSize',10); hold on;
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'CORE_byLME_obs_fished_scatter_all_yrs_P.png'])

%% D
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(Dlme_scatch_all(i,:),Dlme_fcatch_all(i,:),'.','color',[0 0.6 0],'MarkerSize',10); hold on;
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'CORE_byLME_obs_fished_scatter_all_yrs_D.png'])

%% Time series Plots 
% Global
figure(7)
subplot(2,2,1)
plot(years,sum(Flme_scatch_all),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,sum(Flme_fcatch_all),'r','LineWidth',2)
ylabel('FEISTY catch (MT km^-^2)')
title('Forage')

subplot(2,2,2)
plot(years,sum(Plme_scatch_all),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,sum(Plme_fcatch_all),'b','LineWidth',2)
title('Large pelagic')

subplot(2,2,3)
plot(years,sum(Dlme_scatch_all),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,sum(Dlme_fcatch_all),'color',[0 0.6 0],'LineWidth',2)
xlabel('year')
ylabel('FEISTY catch (MT km^-^2)')
title('Demersal')

subplot(2,2,4)
plot(years,sum(slme_scatch_all),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,sum(Alme_fcatch_all),'k','LineWidth',2)
xlabel('year')
title('All')
print('-dpng',[ppath 'CORE_obs_fished_sumlme_ts.png'])

%% log10 Global
figure(8)
subplot(2,2,1)
plot(years,log10(sum(Flme_scatch_all)),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,log10(sum(Flme_fcatch_all)),'r','LineWidth',2)
ylabel('log_1_0 FEISTY catch (MT km^-^2)')
title('Forage')

subplot(2,2,2)
plot(years,log10(sum(Plme_scatch_all)),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,log10(sum(Plme_fcatch_all)),'b','LineWidth',2)
title('Large pelagic')

subplot(2,2,3)
plot(years,log10(sum(Dlme_scatch_all)),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,log10(sum(Dlme_fcatch_all)),'color',[0 0.6 0],'LineWidth',2)
xlabel('year')
ylabel('log_1_0 FEISTY catch (MT km^-^2)')
title('Demersal')

subplot(2,2,4)
plot(years,log10(sum(slme_scatch_all)),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,log10(sum(Alme_fcatch_all)),'k','LineWidth',2)
xlabel('year')
title('All')
print('-dpng',[ppath 'CORE_obs_fished_sumlme_ts_log10.png'])

%% TS by LME
% all 
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(years,slme_scatch_all(i,:),'color',[0.5 0.5 0.5],'LineWidth',1); hold on
        plot(years,Alme_fcatch_all(i,:),'k','LineWidth',1)
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'CORE_byLME_obs_fished_ts_Allfish.png'])

%% F
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(years,Flme_scatch_all(i,:),'color',[0.5 0.5 0.5],'LineWidth',1); hold on
        plot(years,Flme_fcatch_all(i,:),'r','LineWidth',1)
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'CORE_byLME_obs_fished_ts_F.png'])

%% P
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(years,Plme_scatch_all(i,:),'color',[0.5 0.5 0.5],'LineWidth',1); hold on
        plot(years,Plme_fcatch_all(i,:),'b','LineWidth',1)
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'CORE_byLME_obs_fished_ts_P.png'])

%% D
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(years,Dlme_scatch_all(i,:),'color',[0.5 0.5 0.5],'LineWidth',1); hold on
        plot(years,Dlme_fcatch_all(i,:),'color',[0 0.6 0],'LineWidth',1)
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'CORE_byLME_obs_fished_ts_D.png'])


