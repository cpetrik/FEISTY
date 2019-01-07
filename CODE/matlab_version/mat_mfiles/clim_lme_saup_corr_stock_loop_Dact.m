%POEM catch vs. SAUP catch by LME
%Use same methods as Stock et al. 2017 to reduce SAUP dataset
%Diff k values for temp-dep

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/CSV/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/Bio_rates/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'LME_clim_temp.mat']);

%SAUP
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);
l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);
sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

% "lfile" never changes, has lme areas
lfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
lpath = ['/Volumes/GFDL/NC/Matlab_new_size/' lfile '/'];
load([lpath 'LME_clim_fished03_' lfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;
AREA_OCN = max(area,1);

% POEM file info
harv = 'All_fish03';

%% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')
% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

x=-8:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% loop over Demersal metabolism fraction
kays = 0.5:0.1:1;
skays={'.5','.6','.7','.8','.9','1'};

r   = NaN*ones(length(kays),5);
rmse = NaN*ones(length(kays),5);
F   = NaN*ones(length(kays),5);

for i=1:length(kays)
    Dact=kays(i);
    tdact = num2str(1000+int64(100*Dact));
    if (Dact<1)
        cfile = ['Dc_enc70-b200_m4-b175-k08-Dac',tdact(2:end),...
        '_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100'];
        fpath=['/Volumes/GFDL/CSV/Matlab_new_size/' cfile '/'];
    else
        cfile = 'Dc_enc70-b200_m4-b175-k08_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
        fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
    end
    ppath = [pp cfile '/'];
    load([fpath 'LME_clim_fished_',harv,'_' cfile '.mat'],'lme_mcatch');
    
    
    %% POEM LME biomass in MT
    plme_Fmcatch = lme_mcatch(:,1) * 1e-6;
    plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
    plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
    % MT/km2
    plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
    plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
    plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;
    
    l10p=log10(plme_Fmcatch+plme_Pmcatch+plme_Dmcatch);
    l10pF=log10(plme_Fmcatch);
    l10pP=log10(plme_Pmcatch);
    l10pD=log10(plme_Dmcatch);
    
    pFracPD = plme_Pmcatch ./ (plme_Pmcatch + plme_Dmcatch);
    
    %% Stats
    % Drop Arctic, Antarctic, Hawaii, Australia -------------------------
    %r
    r(i,1)=corr(l10s(keep),l10p(keep));
    r(i,2)=corr(l10sF(keep),l10pF(keep));
    r(i,3)=corr(l10sP(keep),l10pP(keep));
    r(i,4)=corr(l10sD(keep),l10pD(keep));
    r(i,5)=corr(sFracPD(keep),pFracPD(keep));
    
    %root mean square error
    o=l10s(keep);
    p=l10p(keep);
    n = length(o);
    num=nansum((p-o).^2);
    rmse(i,1) = sqrt(num/n);
    
    o=l10sF(keep);
    p=l10pF(keep);
    n = length(o);
    num=nansum((p-o).^2);
    rmse(i,2) = sqrt(num/n);
    
    o=l10sP(keep);
    p=l10pP(keep);
    n = length(o);
    num=nansum((p-o).^2);
    rmse(i,3) = sqrt(num/n);
    
    o=l10sD(keep);
    p=l10pD(keep);
    n = length(o);
    num=nansum((p-o).^2);
    rmse(i,4) = sqrt(num/n);
    
    o=sFracPD(keep);
    p=pFracPD(keep);
    n = length(o);
    num=nansum((p-o).^2);
    rmse(i,5) = sqrt(num/n);
    
    %Fmed
    F(i,1)=10^(median(l10s(keep)-l10p(keep)));
    F(i,2)=10^(median(l10sF(keep)-l10pF(keep)));
    F(i,3)=10^(median(l10sP(keep)-l10pP(keep)));
    F(i,4)=10^(median(l10sD(keep)-l10pD(keep)));
    F(i,5)=10^(median(sFracPD(keep)-pFracPD(keep)));
    
    %% Plot corr
    f1=figure(1);
    subplot(2,3,i)
    plot(x,x,'--k');hold on;
    for j=1:length(keep)
        lme=keep(j);
        plot(l10s(lme),l10p(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    end
    text(-2.5,1.7,['r = ' sprintf('%2.2f',r(i,1))])
    text(-2.5,1.3,['RMSE = ' sprintf('%2.2f',rmse(i,1))])
    axis([-3 2 -3 2])
    if (i>6)
        xlabel('SAU')
    end
    ylabel('POEM')
    if (i==2)
        title(['All ' num2str(Dact)])
    else
        title(num2str(Dact))
    end
    stamp('')
    
    f2=figure(2);
    subplot(2,3,i)
    plot(x,x,'--k');hold on;
    for j=1:length(keep)
        lme=keep(j);
        plot(l10sF(lme),l10pF(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    end
    text(-5.5,1.7,['r = ' sprintf('%2.2f',r(i,2))])
    text(-5.5,1.1,['RMSE = ' sprintf('%2.2f',rmse(i,2))])
    axis([-6 2 -6 2])
    if (i>6)
        xlabel('SAU')
    end
    ylabel('POEM')
    if (i==2)
        title(['Forage ' num2str(Dact)])
    else
        title(num2str(Dact))
    end
    stamp('')
    
    f3=figure(3);
    subplot(2,3,i)
    plot(x,x,'--k');hold on;
    for j=1:length(keep)
        lme=keep(j);
        plot(l10sP(lme),l10pP(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    end
    %plot(l10sP,l10pP,'.k','MarkerSize',10); hold on;
    text(-7.5,1.5,['r = ' sprintf('%2.2f',r(i,3))])
    text(-7.5,0.5,['RMSE = ' sprintf('%2.2f',rmse(i,3))])
    axis([-8 2 -8 2])
    if (i>6)
        xlabel('SAU')
    end
    ylabel('POEM')
    if (i==2)
        title(['Large Pelagics ' num2str(Dact)])
    else
        title(num2str(Dact))
    end
    stamp('')

    f4=figure(4);
    subplot(2,3,i)
    plot(x,x,'--k');hold on;
    for j=1:length(keep)
        lme=keep(j);
        plot(l10sD(lme),l10pD(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    end
    %plot(l10sD,l10pD,'.k','MarkerSize',10); hold on;
    text(-1.75,1.7,['r = ' sprintf('%2.2f',r(i,4))])
    text(-1.75,1.3,['RMSE = ' sprintf('%2.2f',rmse(i,4))])
    axis([-2 2 -2 2])
    if (i>6)
        xlabel('SAU')
    end
    ylabel('POEM')
    if (i==2)
        title(['Demersals ' num2str(Dact)])
    else
        title(num2str(Dact))
    end
    stamp('')
    
    f5=figure(5);
    subplot(2,3,i)
    plot(x,x,'--k');hold on;
    for j=1:length(keep)
        lme=keep(j);
        plot(sFracPD(lme),pFracPD(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    end
    text(0.7,0.35,['r = ' sprintf('%2.2f',r(i,5))])
    text(0.7,0.25,['RMSE = ' sprintf('%2.2f',rmse(i,5))])
    axis([0 1 0 1])
    if (i>6)
        xlabel('SAU')
    end
    ylabel('POEM')
    if (i==2)
        title(['P/(P+D) ' num2str(Dact)])
    else
        title(num2str(Dact))
    end
    stamp('')

end
print(f1,'-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_Dact_tests_scatterAll.png'])
print(f2,'-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_Dact_tests_scatterF.png'])
print(f3,'-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_Dact_tests_scatterP.png'])
print(f4,'-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_Dact_tests_scatterD.png'])
print(f5,'-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_Dact_tests_scatterPD.png'])
    
%% Plots

figure(6)
subplot(2,2,1)
plot(kays,r(:,1),'k','LineWidth',2); hold on;
plot(kays,r(:,2),'color',[0.97647 0.45098 0.023529],'LineWidth',2)
plot(kays,r(:,3),'b','LineWidth',2); hold on;
plot(kays,r(:,4),'color',[0.082353 0.6902 0.10196],'LineWidth',2)
plot(kays,r(:,5),'color',[0.49412 0.11765 0.61176],'LineWidth',2)
xlabel('Dact')
ylabel('r')
xlim([0.5 1])

subplot(2,2,2)
plot(kays,rmse(:,1),'k','LineWidth',2); hold on;
plot(kays,rmse(:,2),'color',[0.97647 0.45098 0.023529],'LineWidth',2)
plot(kays,rmse(:,3),'b','LineWidth',2); hold on;
plot(kays,rmse(:,4),'color',[0.082353 0.6902 0.10196],'LineWidth',2)
plot(kays,rmse(:,5),'color',[0.49412 0.11765 0.61176],'LineWidth',2)
xlabel('Dact')
ylabel('rmse')
xlim([0.5 1])
legend('All','F','P','D','P:D')
legend('location','northwest')

subplot(2,2,3)
plot(kays,F(:,1),'k','LineWidth',2); hold on;
plot(kays,F(:,2),'color',[0.97647 0.45098 0.023529],'LineWidth',2)
plot(kays,F(:,3),'b','LineWidth',2); hold on;
plot(kays,F(:,4),'color',[0.082353 0.6902 0.10196],'LineWidth',2)
plot(kays,F(:,5),'color',[0.49412 0.11765 0.61176],'LineWidth',2)
xlabel('Dact')
ylabel('Fmed')
xlim([0.5 1])
stamp('')
print('-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_Dact_tests_stats_line.png'])

%%
figure(7)
subplot(2,2,1)
plot(kays,r(:,5),'k','LineWidth',2); hold on;
xlabel('Dact')
ylabel('r')
xlim([0.5 1])
title('P/(P+D)')

subplot(2,2,2)
plot(kays,rmse(:,5),'k','LineWidth',2); hold on;
xlabel('Dact')
ylabel('rmse')
xlim([0.5 1])

subplot(2,2,3)
plot(kays,F(:,5),'k','LineWidth',2); hold on;
xlabel('Dact')
ylabel('Fmed')
xlim([0.5 1])
stamp('')
print('-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_Dact_tests_PDstats_line.png'])

