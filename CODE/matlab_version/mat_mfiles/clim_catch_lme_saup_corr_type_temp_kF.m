%POEM catch vs. SAUP catch by LME

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

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

lfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
lpath = [dp lfile '/'];
load([lpath 'LME_clim_fished03_' lfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

% plot info
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

%%

Flme_catch_all = nansum(Flme_wcatch,3);
Plme_catch_all = nansum(Plme_wcatch,3);
Dlme_catch_all = nansum(Dlme_wcatch,3);

%1956-2005 SAUP average
id = find(yr>1955 & yr<=2005);

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
%Top 10 yrs SAUP
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

%% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);

keep = [1:9,11:34,36:38,46:54,59:60,62,65];

x=-6:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% check with a plot
figure(100)
for i=1:66
    plot(i,i,'.','color',tmap(tid(i,2),:)); hold on;
end

%% Loop over POEM params
bpow = 0.175;
kays = 0.0405:0.01:0.125;
Fish = 0.1:0.1:0.6;

for b=1:length(Fish)
    for k = 1:6;%length(kays)
        frate = Fish(b); 
        kt = kays(k);
        tkfn = num2str(100+int64(100*kt));
        tbfn = num2str(100+int64(100*bpow));
        tfish = num2str(100+int64(10*frate));
        
        cfile = ['Dc_enc70_cmax-metab20_b',tbfn(2:end),'_k',tkfn(2:end),'_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100'];
        
        charv = ['All_fish',tfish(2:end)];
        harv = tfish(2:end);
        
        dpath = [dp cfile '/'];
        ppath = [pp cfile '/'];
        
        load([dpath 'LME_clim_fished',harv,'_' cfile '.mat']);
        
        close all
        
        %% POEM LME biomass in MT
        plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
        plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
        plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
        plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
        % MT/km2
        plme_mcatch = plme_mcatch ./ lme_area_km2;
        plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
        plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
        plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;
        
        %% Plots
        l10s=log10(slme_mcatch10+eps);
        l10p=log10(plme_mcatch);
        l10sF=log10(Flme_mcatch10+eps);
        l10pF=log10(plme_Fmcatch);
        l10sP=log10(Plme_mcatch10+eps);
        l10pP=log10(plme_Pmcatch);
        l10sD=log10(Dlme_mcatch10+eps);
        l10pD=log10(plme_Dmcatch);
        
        figure(1)
        subplot(2,2,1)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,'--b'); hold on;
        plot(x,x2l,'--b'); hold on;
        plot(x,x5h,'--r'); hold on;
        plot(x,x5l,'--r'); hold on;
        for i=1:66
            plot(l10s(i),l10p(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
        end
        axis([-6 2 -6 2])
        xlabel('SAUP mean of top 10 years')
        ylabel('POEM mean of Climatology')
        title('Mean catch')
        
        subplot(2,2,2)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,'--b'); hold on;
        plot(x,x2l,'--b'); hold on;
        plot(x,x5h,'--r'); hold on;
        plot(x,x5l,'--r'); hold on;
        for i=1:66
            plot(l10sF(i),l10pF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
        end
        axis([-6 2 -6 2])
        xlabel('SAUP F catch (log10 MT km^-^2)')
        ylabel('POEM F catch (log10 MT km^-^2)')
        title('Mean F catch')
        
        subplot(2,2,3)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,'--b'); hold on;
        plot(x,x2l,'--b'); hold on;
        plot(x,x5h,'--r'); hold on;
        plot(x,x5l,'--r'); hold on;
        for i=1:66
            plot(l10sP(i),l10pP(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
        end
        axis([-6 2 -6 2])
        xlabel('SAUP P catch (log10 MT km^-^2)')
        ylabel('POEM P catch (log10 MT km^-^2)')
        title('Mean P catch')
        
        subplot(2,2,4)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,'--b'); hold on;
        plot(x,x2l,'--b'); hold on;
        plot(x,x5h,'--r'); hold on;
        plot(x,x5l,'--r'); hold on;
        for i=1:66
            plot(l10sD(i),l10pD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
        end
        axis([-6 2 -6 2])
        xlabel('SAUP D catch (log10 MT km^-^2)')
        ylabel('POEM D catch (log10 MT km^-^2)')
        title('Mean D catch')
        stamp(cfile)
        print('-dpng',[ppath 'Clim_fished',harv,'_SAUP_log10catch_comp_types_temp.png'])
        
        
        
        %% Drop Arctic, Antarctic, Hawaii, Australia -------------------------
        % Stats
        %r
        rall=corr(l10s(keep),l10p(keep));
        rF=corr(l10sF(keep),l10pF(keep));
        rP=corr(l10sP(keep),l10pP(keep));
        rD=corr(l10sD(keep),l10pD(keep));
        
        %root mean square error
        o=l10s(keep);
        p=l10p(keep);
        n = length(o);
        num=nansum((p-o).^2);
        rmse = sqrt(num/n);
        
        o=l10sF(keep);
        p=l10pF(keep);
        n = length(o);
        num=nansum((p-o).^2);
        rmseF = sqrt(num/n);
        
        o=l10sP(keep);
        p=l10pP(keep);
        n = length(o);
        num=nansum((p-o).^2);
        rmseP = sqrt(num/n);
        
        o=l10sD(keep);
        p=l10pD(keep);
        n = length(o);
        num=nansum((p-o).^2);
        rmseD = sqrt(num/n);
        
        
        figure(11)
        subplot(2,2,1)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,'--b'); hold on;
        plot(x,x2l,'--b'); hold on;
        plot(x,x5h,'--r'); hold on;
        plot(x,x5l,'--r'); hold on;
        for i=1:length(keep)
            lme=keep(i);
            plot(l10s(lme),l10p(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
        end
        text(-3.5,1.5,['r = ' num2str(rall)])
        text(-3.5,1.0,['RMSE = ' num2str(rmse)])
        axis([-4 2 -4 2])
        xlabel('SAUP catch (log10 MT km^-^2) mean of top 10 years')
        ylabel('Total mean of Climatology')
        title('Mean catch without Polar and Australia')
        
        subplot(2,2,2)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,'--b'); hold on;
        plot(x,x2l,'--b'); hold on;
        plot(x,x5h,'--r'); hold on;
        plot(x,x5l,'--r'); hold on;
        for i=1:length(keep)
            lme=keep(i);
            plot(l10sF(lme),l10pF(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
        end
        text(-5.5,1.5,['r = ' num2str(rF)])
        text(-5.5,1.0,['RMSE = ' num2str(rmseF)])
        axis([-6 2 -6 2])
        xlabel('SAUP F catch (log10 MT km^-^2)')
        ylabel('POEM F catch (log10 MT km^-^2)')
        title('Mean F catch')
        
        subplot(2,2,3)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,'--b'); hold on;
        plot(x,x2l,'--b'); hold on;
        plot(x,x5h,'--r'); hold on;
        plot(x,x5l,'--r'); hold on;
        for i=1:length(keep)
            lme=keep(i);
            plot(l10sP(lme),l10pP(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
        end
        text(-5.5,1.5,['r = ' num2str(rP)])
        text(-5.5,1.0,['RMSE = ' num2str(rmseP)])
        axis([-6 2 -6 2])
        xlabel('SAUP P catch (log10 MT km^-^2)')
        ylabel('POEM P catch (log10 MT km^-^2)')
        title('Mean P catch')
        
        subplot(2,2,4)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,'--b'); hold on;
        plot(x,x2l,'--b'); hold on;
        plot(x,x5h,'--r'); hold on;
        plot(x,x5l,'--r'); hold on;
        for i=1:length(keep)
            lme=keep(i);
            plot(l10sD(lme),l10pD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
        end
        text(-5.5,1.5,['r = ' num2str(rD)])
        text(-5.5,1.0,['RMSE = ' num2str(rmseD)])
        axis([-6 2 -6 2])
        xlabel('SAUP D catch (log10 MT km^-^2)')
        ylabel('POEM D catch (log10 MT km^-^2)')
        title('Mean D catch')
        stamp(cfile)
        print('-dpng',[ppath 'Clim_fished',harv,'_SAUP_log10catch_comp_types_LELC_temp.png'])
        
        
        
    end
end


