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

%use weighted catches
load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
    'Flme_wcatch','Dlme_wcatch','Plme_wcatch');

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')

AREA_OCN = max(area,1);

% "lfile" never changes, has lme areas
lfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
lpath = ['/Volumes/GFDL/NC/Matlab_new_size/' lfile '/'];
load([lpath 'LME_clim_fished03_' lfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

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

% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

x=-8:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% SAUP
Flme_catch_all = nansum(Flme_wcatch,3);
Plme_catch_all = nansum(Plme_wcatch,3);
Dlme_catch_all = nansum(Dlme_wcatch,3);

%1950-2006 SAUP average
id = find(yr>1950 & yr<=2006);

slme_mcatch = nanmean(lme_catch(id,:));
slme_mcatch = slme_mcatch';
Fslme_mcatch = nanmean(Flme_catch_all(id,:));
Fslme_mcatch = Fslme_mcatch';
Pslme_mcatch = nanmean(Plme_catch_all(id,:));
Pslme_mcatch = Pslme_mcatch';
Dslme_mcatch = nanmean(Dlme_catch_all(id,:));
Dslme_mcatch = Dslme_mcatch';

slme_mcatch10 = NaN*ones(size(slme_mcatch));
Flme_mcatch10 = NaN*ones(size(slme_mcatch));
Plme_mcatch10 = NaN*ones(size(slme_mcatch));
Dlme_mcatch10 = NaN*ones(size(slme_mcatch));

%Top 10 yrs by LME SAUP
for i=1:66
    [sort_lme_catch,ix] = sort(lme_catch(:,i),'descend');
    sort_Flme_catch = Flme_catch_all(ix,i);
    sort_Plme_catch = Plme_catch_all(ix,i);
    sort_Dlme_catch = Dlme_catch_all(ix,i);
    slme_mcatch10(i) = nanmean(sort_lme_catch(1:10));
    Flme_mcatch10(i) = nanmean(sort_Flme_catch(1:10));
    Plme_mcatch10(i) = nanmean(sort_Plme_catch(1:10));
    Dlme_mcatch10(i) = nanmean(sort_Dlme_catch(1:10));
end

% MT/km2
slme_mcatch10 = slme_mcatch10 ./ lme_area_km2;
Flme_mcatch10 = Flme_mcatch10 ./ lme_area_km2;
Plme_mcatch10 = Plme_mcatch10 ./ lme_area_km2;
Dlme_mcatch10 = Dlme_mcatch10 ./ lme_area_km2;

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

%% loop over temp scaling
kays = 0.0405:0.01:0.1205;
skays={'.04','.05','.06','.07','.08','.09','.10','.11','.12'};

r   = NaN*ones(length(kays),5);
rmse = NaN*ones(length(kays),5);
F   = NaN*ones(length(kays),5);

for i=1:length(kays)
    kt=kays(i);
    tkfn = num2str(100+int64(100*kt));
    cfile = ['Dc_enc70-b200_cm20_m-b175-k',tkfn(2:end),...
        '_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100'];
    fpath=['/Volumes/GFDL/CSV/Matlab_new_size/' cfile '/'];
    ppath = [pp cfile '/'];
    load([fpath 'LME_clim_',harv,'_loop.mat'],'lme_mcatch');
    
    
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
    subplot(3,3,i)
    plot(x,x,'--k');hold on;
    for j=1:length(keep)
        lme=keep(j);
        plot(l10s(lme),l10p(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    end
    text(-1.75,1.7,['r = ' sprintf('%2.2f',r(i,1))])
    text(-1.75,1.3,['RMSE = ' sprintf('%2.2f',rmse(i,1))])
    axis([-3 2 -3 2])
    if (i>6)
        xlabel('SAU')
    end
    ylabel('POEM')
    if (i==2)
        title('All')
    else
        title(['k=0' skays{i}])
    end
    stamp('')
    
    f2=figure(2);
    subplot(3,3,i)
    plot(x,x,'--k');hold on;
    for j=1:length(keep)
        lme=keep(j);
        plot(l10sF(lme),l10pF(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    end
    text(-7.5,1.5,['r = ' sprintf('%2.2f',r(i,2))])
    text(-7.5,0.5,['RMSE = ' sprintf('%2.2f',rmse(i,2))])
    axis([-8 2 -8 2])
    if (i>6)
        xlabel('SAU')
    end
    ylabel('POEM')
    if (i==2)
        title('Forage')
    else
        title(['k=0' skays{i}])
    end
    stamp('')
    
    f3=figure(3);
    subplot(3,3,i)
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
        title('Large Pelagics')
    else
        title(['k=0' skays{i}])
    end
    stamp('')

    f4=figure(4);
    subplot(3,3,i)
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
        title('Demersals')
    else
        title(['k=0' skays{i}])
    end
    stamp('')
    
    f5=figure(5);
    subplot(3,3,i)
    plot(x,x,'--k');hold on;
    for j=1:length(keep)
        lme=keep(j);
        plot(sFracPD(lme),pFracPD(lme),'.','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    end
    text(0.8,0.35,['r = ' sprintf('%2.2f',r(i,5))])
    text(0.8,0.25,['RMSE = ' sprintf('%2.2f',rmse(i,5))])
    axis([0 1 0 1])
    if (i>6)
        xlabel('SAU')
    end
    ylabel('POEM')
    if (i==2)
        title('P/(P+D)')
    else
        title(['k=0' skays{i}])
    end
    stamp('')

end
print(f1,'-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_kt_tests_scatterAll.png'])
print(f2,'-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_kt_tests_scatterF.png'])
print(f3,'-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_kt_tests_scatterP.png'])
print(f4,'-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_kt_tests_scatterD.png'])
print(f5,'-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_kt_tests_scatterPD.png'])
    
%% Plots

figure(6)
subplot(2,2,1)
plot(kays,r(:,1),'k','LineWidth',2); hold on;
plot(kays,r(:,2),'color',[0.97647 0.45098 0.023529],'LineWidth',2)
plot(kays,r(:,3),'b','LineWidth',2); hold on;
plot(kays,r(:,4),'color',[0.082353 0.6902 0.10196],'LineWidth',2)
plot(kays,r(:,5),'color',[0.49412 0.11765 0.61176],'LineWidth',2)
xlabel('k')
ylabel('r')
xlim([0.04 0.12])

subplot(2,2,2)
plot(kays,rmse(:,1),'k','LineWidth',2); hold on;
plot(kays,rmse(:,2),'color',[0.97647 0.45098 0.023529],'LineWidth',2)
plot(kays,rmse(:,3),'b','LineWidth',2); hold on;
plot(kays,rmse(:,4),'color',[0.082353 0.6902 0.10196],'LineWidth',2)
plot(kays,rmse(:,5),'color',[0.49412 0.11765 0.61176],'LineWidth',2)
xlabel('k')
ylabel('rmse')
xlim([0.04 0.12])
legend('All','F','P','D','P:D')
legend('location','northwest')

subplot(2,2,3)
plot(kays,F(:,1),'k','LineWidth',2); hold on;
plot(kays,F(:,2),'color',[0.97647 0.45098 0.023529],'LineWidth',2)
plot(kays,F(:,3),'b','LineWidth',2); hold on;
plot(kays,F(:,4),'color',[0.082353 0.6902 0.10196],'LineWidth',2)
plot(kays,F(:,5),'color',[0.49412 0.11765 0.61176],'LineWidth',2)
xlabel('k')
ylabel('Fmed')
xlim([0.04 0.12])
stamp('')
print('-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_kt_tests_stats_line.png'])

%%
figure(7)
subplot(2,2,1)
plot(kays,r(:,5),'k','LineWidth',2); hold on;
xlabel('k')
ylabel('r')
xlim([0.04 0.12])
title('P/(P+D)')

subplot(2,2,2)
plot(kays,rmse(:,5),'k','LineWidth',2); hold on;
xlabel('k')
ylabel('rmse')
xlim([0.04 0.12])

subplot(2,2,3)
plot(kays,F(:,5),'k','LineWidth',2); hold on;
xlabel('k')
ylabel('Fmed')
xlim([0.04 0.12])
stamp('')
print('-dpng',[pp 'Clim_',harv,'_SAUP_comp_Stock_LELC_kt_tests_PDstats_line.png'])

