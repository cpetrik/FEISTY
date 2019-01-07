%POEM catch vs. SAUP catch by LME

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

load([cpath 'hindcast_gridspec.mat'],'AREA_OCN','geolat_t','geolon_t');
load([cpath 'lme_mask_esm2m.mat']);
load([cpath 'LME_clim_temp.mat']);
grid = csvread([cpath 'grid_csv.csv']);

%use weighted catches
load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
    'Flme_wcatch','Dlme_wcatch','Plme_wcatch');

%area of LMEs
lfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
lpath = [dp lfile '/'];
load([lpath 'LME_clim_fished03_' lfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')
cmap1(1,:)=[1 1 1];
cmap1(2,:)=cmap_ppt(1,:);
cmap1(3,:)=cmap_ppt(3,:);
cmap1(4,:)=cmap_ppt(5,:);

cmap2(1,:)=cmap_ppt(1,:);
cmap2(2,:)=cmap_ppt(3,:);
cmap2(3,:)=cmap_ppt(5,:);

cmap3(1,:)=[0 0 0];
cmap3(2,:)=cmap_ppt(1,:);
cmap3(3,:)=cmap_ppt(3,:);
cmap3(4,:)=cmap_ppt(5,:);

cmap4(1,:)=cmap_ppt(3,:);
cmap4(2,:)=cmap_ppt(1,:);
cmap4(3,:)=cmap_ppt(5,:);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

% plot info
[ni,nj]=size(geolon_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];
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

%% Loop over POEM params
harv = '03';

% kays = 0.0405:0.01:0.0916;
% bees = 0.125:0.025:0.25; %bees = 0.1:0.05:0.35;
% 
% rall = NaN*ones(length(kays),length(bees));
% rF = rall;
% rP = rall;
% rD = rall;
% rmse = rall;
% rmseF = rall;
% rmseP = rall;
% rmseD = rall;
% 
% 
% for k = 1:length(kays)
%     kt = kays(k);
%     
%     for b = 1:length(bees)
%         bpow = bees(b);
%         
%         tkfn = num2str(100+int64(100*kt));
%         tbfn = num2str(100+int64(100*bpow));
%         
%         cfile = ['Dc_enc70_cmax-metab20_b',tbfn(2:end),'_k',tkfn(2:end),'_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100'];
%         
        cfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC100_lgRE00100_mdRE00100';
        dpath = [dp cfile '/'];
        ppath = [pp cfile '/'];
        
        load([dpath 'LME_hindcast_fished',harv,'_' cfile '.mat']);
        
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
        
        %log10 Difference
        diff_catch = (log10(plme_mcatch) - log10(slme_mcatch10)) ./ log10((slme_mcatch10+eps));
        Fdiff_catch = (log10(plme_Fmcatch) - log10(Flme_mcatch10)) ./ log10((Flme_mcatch10+eps));
        Pdiff_catch = (log10(plme_Pmcatch) - log10(Plme_mcatch10)) ./ log10((Plme_mcatch10+eps));
        Ddiff_catch = (log10(plme_Dmcatch) - log10(Dlme_mcatch10)) ./ log10((Dlme_mcatch10+eps));
        
        code = [1:66]';
        T = table(code,slme_mcatch10,plme_mcatch,Flme_mcatch10,plme_Fmcatch,...
            Plme_mcatch10,plme_Pmcatch,Dlme_mcatch10,plme_Dmcatch,'VariableNames',...
            {'lme','saupAll','poemAll','saupF','poemF','saupP','poemP','saupD','poemD'});
        writetable(T,[dpath 'LME_saup_catch_spinup_fished',harv,'.csv']);
        
        %% Plots
        l10s=log10(slme_mcatch10+eps);
        l10p=log10(plme_mcatch);
        l10sF=log10(Flme_mcatch10+eps);
        l10pF=log10(plme_Fmcatch);
        l10sP=log10(Plme_mcatch10+eps);
        l10pP=log10(plme_Pmcatch);
        l10sD=log10(Dlme_mcatch10+eps);
        l10pD=log10(plme_Dmcatch);
        
        figure(10)
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
        xlabel('SAUP catch mean of top 10 years')
        ylabel('POEM catch mean of Historic 1996-2005')
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
        print('-dpng',[ppath 'Hindcast_fished',harv,'_SAUP_log10catch_comp_types_temp.png'])
        
        
        %% Plot by region
        
        nwa = [6:9,18];
        nea = [19:26,59:60];
        nep = [1:4,10:11,54:55,64:65];
        nwp = [47:53,56];
        sam = 13:17;
        afr = [27:31,33];
        aus = 39:46;
        ind = [32,34:38];
        crb = [5,12];
        
%         %F
%         figure(5)
%         subplot(2,2,1)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Flme_mcatch10(nwa)),log10(plme_Fmcatch(nwa)),'.c','MarkerSize',25); hold on;
%         axis([0 8 0 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('F NW Atl 1996-2005 mean catch')
%         for i=1:length(nwa)
%             text(log10(Flme_mcatch10(nwa(i))),log10(plme_Fmcatch(nwa(i))),num2str(nwa(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,2,2)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Flme_mcatch10(nea)),log10(plme_Fmcatch(nea)),'.c','MarkerSize',25); hold on;
%         axis([0 8 0 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('F NE Atl')
%         for i=1:length(nea)
%             text(log10(Flme_mcatch10(nea(i))),log10(plme_Fmcatch(nea(i))),num2str(nea(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,2,3)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Flme_mcatch10(nwp)),log10(plme_Fmcatch(nwp)),'.c','MarkerSize',25); hold on;
%         axis([0 8 0 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('F NW Pac')
%         for i=1:length(nwp)
%             text(log10(Flme_mcatch10(nwp(i))),log10(plme_Fmcatch(nwp(i))),num2str(nwp(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,2,4)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Flme_mcatch10(nep)),log10(plme_Fmcatch(nep)),'.c','MarkerSize',25); hold on;
%         axis([0 8 0 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('F NE Pac')
%         for i=1:length(nep)
%             text(log10(Flme_mcatch10(nep(i))),log10(plme_Fmcatch(nep(i))),num2str(nep(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         print('-dpng',[ppath 'Historic_fished',harv,'_SAUP10_log10catch_compF_AtlPac.png'])
%         
%         %% P
%         figure(6)
%         subplot(2,2,1)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Plme_mcatch10(nwa)),log10(plme_Pmcatch(nwa)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('P NW Atl 1996-2005 mean catch')
%         for i=1:length(nwa)
%             text(log10(Plme_mcatch10(nwa(i))),log10(plme_Pmcatch(nwa(i))),num2str(nwa(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,2,2)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Plme_mcatch10(nea)),log10(plme_Pmcatch(nea)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('P NE Atl')
%         for i=1:length(nea)
%             text(log10(Plme_mcatch10(nea(i))),log10(plme_Pmcatch(nea(i))),num2str(nea(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,2,3)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Plme_mcatch10(nwp)),log10(plme_Pmcatch(nwp)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('P NW Pac')
%         for i=1:length(nwp)
%             text(log10(Plme_mcatch10(nwp(i))),log10(plme_Pmcatch(nwp(i))),num2str(nwp(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,2,4)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Plme_mcatch10(nep)),log10(plme_Pmcatch(nep)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('P NE Pac')
%         for i=1:length(nep)
%             text(log10(Plme_mcatch10(nep(i))),log10(plme_Pmcatch(nep(i))),num2str(nep(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         print('-dpng',[ppath 'Historic_fished',harv,'_SAUP10_log10catch_compP_AtlPac.png'])
%         
%         %D
%         figure(7)
%         subplot(2,2,1)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Dlme_mcatch10(nwa)),log10(plme_Dmcatch(nwa)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('D NW Atl 1996-2005 mean catch')
%         for i=1:length(nwa)
%             text(log10(Dlme_mcatch10(nwa(i))),log10(plme_Dmcatch(nwa(i))),num2str(nwa(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,2,2)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Dlme_mcatch10(nea)),log10(plme_Dmcatch(nea)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('D NE Atl')
%         for i=1:length(nea)
%             text(log10(Dlme_mcatch10(nea(i))),log10(plme_Dmcatch(nea(i))),num2str(nea(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,2,3)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Dlme_mcatch10(nwp)),log10(plme_Dmcatch(nwp)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('D NW Pac')
%         for i=1:length(nwp)
%             text(log10(Dlme_mcatch10(nwp(i))),log10(plme_Dmcatch(nwp(i))),num2str(nwp(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,2,4)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Dlme_mcatch10(nep)),log10(plme_Dmcatch(nep)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('D NE Pac')
%         for i=1:length(nep)
%             text(log10(Dlme_mcatch10(nep(i))),log10(plme_Dmcatch(nep(i))),num2str(nep(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         print('-dpng',[ppath 'Historic_fished',harv,'_SAUP10_log10catch_compD_AtlPac.png'])
%         
%         %%
%         figure(8)
%         subplot(2,3,1)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Flme_mcatch10(afr)),log10(plme_Fmcatch(afr)),'.c','MarkerSize',25); hold on;
%         axis([0 8 0 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('F Africa 1996-2005 mean catch')
%         for i=1:length(afr)
%             text(log10(Flme_mcatch10(afr(i))),log10(plme_Fmcatch(afr(i))),num2str(afr(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,3,2)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Flme_mcatch10(aus)),log10(plme_Fmcatch(aus)),'.c','MarkerSize',25); hold on;
%         axis([0 8 0 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('F Australia')
%         for i=1:length(aus)
%             text(log10(Flme_mcatch10(aus(i))),log10(plme_Fmcatch(aus(i))),num2str(aus(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,3,3)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Flme_mcatch10(ind)),log10(plme_Fmcatch(ind)),'.c','MarkerSize',25); hold on;
%         axis([0 8 0 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('F India-Indonesia')
%         for i=1:length(ind)
%             text(log10(Flme_mcatch10(ind(i))),log10(plme_Fmcatch(ind(i))),num2str(ind(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,3,4)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Flme_mcatch10(crb)),log10(plme_Fmcatch(crb)),'.c','MarkerSize',25); hold on;
%         axis([0 8 0 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('F Caribbean')
%         for i=1:length(crb)
%             text(log10(Flme_mcatch10(crb(i))),log10(plme_Fmcatch(crb(i))),num2str(crb(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,3,5)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Flme_mcatch10(sam)),log10(plme_Fmcatch(sam)),'.c','MarkerSize',25); hold on;
%         axis([0 8 0 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('F S Amer')
%         for i=1:length(sam)
%             text(log10(Flme_mcatch10(sam(i))),log10(plme_Fmcatch(sam(i))),num2str(sam(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         print('-dpng',[ppath 'Historic_fished',harv,'_SAUP10_log10catch_compF_Other.png'])
%         
%         
%         figure(9)
%         subplot(2,3,1)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Plme_mcatch10(afr)),log10(plme_Pmcatch(afr)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('P Africa 1996-2005 mean catch')
%         for i=1:length(afr)
%             text(log10(Plme_mcatch10(afr(i))),log10(plme_Pmcatch(afr(i))),num2str(afr(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,3,2)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Plme_mcatch10(aus)),log10(plme_Pmcatch(aus)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('P Australia')
%         for i=1:length(aus)
%             text(log10(Plme_mcatch10(aus(i))),log10(plme_Pmcatch(aus(i))),num2str(aus(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,3,3)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Plme_mcatch10(ind)),log10(plme_Pmcatch(ind)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('P India-Indonesia')
%         for i=1:length(ind)
%             text(log10(Plme_mcatch10(ind(i))),log10(plme_Pmcatch(ind(i))),num2str(ind(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,3,4)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Plme_mcatch10(crb)),log10(plme_Pmcatch(crb)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('P Caribbean')
%         for i=1:length(crb)
%             text(log10(Plme_mcatch10(crb(i))),log10(plme_Pmcatch(crb(i))),num2str(crb(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,3,5)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Plme_mcatch10(sam)),log10(plme_Pmcatch(sam)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('P S Amer')
%         for i=1:length(sam)
%             text(log10(Plme_mcatch10(sam(i))),log10(plme_Pmcatch(sam(i))),num2str(sam(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         print('-dpng',[ppath 'Historic_fished',harv,'_SAUP10_log10catch_compP_Other.png'])
%         
%         
%         figure(10)
%         subplot(2,3,1)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Dlme_mcatch10(afr)),log10(plme_Dmcatch(afr)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('D Africa 1996-2005 mean catch')
%         for i=1:length(afr)
%             text(log10(Dlme_mcatch10(afr(i))),log10(plme_Dmcatch(afr(i))),num2str(afr(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,3,2)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Dlme_mcatch10(aus)),log10(plme_Dmcatch(aus)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('D Australia')
%         for i=1:length(aus)
%             text(log10(Dlme_mcatch10(aus(i))),log10(plme_Dmcatch(aus(i))),num2str(aus(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,3,3)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Dlme_mcatch10(ind)),log10(plme_Dmcatch(ind)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('D India-Indonesia')
%         for i=1:length(ind)
%             text(log10(Dlme_mcatch10(ind(i))),log10(plme_Dmcatch(ind(i))),num2str(ind(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,3,4)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Dlme_mcatch10(crb)),log10(plme_Dmcatch(crb)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('D Caribbean')
%         for i=1:length(crb)
%             text(log10(Dlme_mcatch10(crb(i))),log10(plme_Dmcatch(crb(i))),num2str(crb(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         
%         subplot(2,3,5)
%         plot(x,x,'--k'); hold on;
%         plot(x,x2h,'--b'); hold on;
%         plot(x,x2l,'--b'); hold on;
%         plot(x,x5h,'--r'); hold on;
%         plot(x,x5l,'--r'); hold on;
%         plot(log10(Dlme_mcatch10(sam)),log10(plme_Dmcatch(sam)),'.c','MarkerSize',25); hold on;
%         axis([2 8 2 8])
%         xlabel('log10 SAUP catch (MT) top 10 yrs')
%         ylabel('log10 POEM catch (MT)')
%         title('D S Amer')
%         for i=1:length(sam)
%             text(log10(Dlme_mcatch10(sam(i))),log10(plme_Dmcatch(sam(i))),num2str(sam(i)),...
%                 'Color','k','HorizontalAlignment','center')
%         end
%         print('-dpng',[ppath 'Historic_fished',harv,'_SAUP10_log10catch_compD_Other.png'])
%         
        %% Drop Arctic, Antarctic, Hawaii, Australia -------------------------
        
        keep = [1:9,11:34,36:38,46:54,59:60,62,65];
        
        % Stats
        %r
        rall(k,b)=corr(l10s(keep),l10p(keep));
        rF(k,b)=corr(l10sF(keep),l10pF(keep));
        rP(k,b)=corr(l10sP(keep),l10pP(keep));
        rD(k,b)=corr(l10sD(keep),l10pD(keep));
        
        %root mean square error
        o=l10s(keep);
        p=l10p(keep);
        n = length(o);
        num=nansum((p-o).^2);
        rmse(k,b) = sqrt(num/n);
        
        o=l10sF(keep);
        p=l10pF(keep);
        n = length(o);
        num=nansum((p-o).^2);
        rmseF(k,b) = sqrt(num/n);
        
        o=l10sP(keep);
        p=l10pP(keep);
        n = length(o);
        num=nansum((p-o).^2);
        rmseP(k,b) = sqrt(num/n);
        
        o=l10sD(keep);
        p=l10pD(keep);
        n = length(o);
        num=nansum((p-o).^2);
        rmseD(k,b) = sqrt(num/n);
        
        figure(11)
        subplot(2,2,1)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,'--b'); hold on;
        plot(x,x2l,'--b'); hold on;
        plot(x,x5h,'--r'); hold on;
        plot(x,x5l,'--r'); hold on;
        plot(l10s(keep),l10p(keep),'.k','MarkerSize',25); hold on;
        text(3.5,7.5,['r = ' num2str(rall(k,b))])
        text(3.5,7.0,['RMSE = ' num2str(rmse(k,b))])
        axis([3 8 3 8])
        xlabel('SAUP catch (MT) mean of top 10 years')
        ylabel('Total mean of Historic 1996-2005')
        title('Mean catch without Polar and Australia')
        
        subplot(2,2,2)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,'--b'); hold on;
        plot(x,x2l,'--b'); hold on;
        plot(x,x5h,'--r'); hold on;
        plot(x,x5l,'--r'); hold on;
        plot(l10sF(keep),l10pF(keep),'.k','MarkerSize',25); hold on;
        text(1,7.5,['r = ' num2str(rF(k,b))])
        text(1,7.0,['RMSE = ' num2str(rmseF(k,b))])
        axis([0 8 0 8])
        xlabel('SAUP F catch (MT)')
        ylabel('POEM F catch (MT)')
        title('Mean F catch')
        
        subplot(2,2,3)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,'--b'); hold on;
        plot(x,x2l,'--b'); hold on;
        plot(x,x5h,'--r'); hold on;
        plot(x,x5l,'--r'); hold on;
        plot(l10sP(keep),l10pP(keep),'.k','MarkerSize',25); hold on;
        text(2,7.5,['r = ' num2str(rP(k,b))])
        text(2,7.0,['RMSE = ' num2str(rmseP(k,b))])
        axis([1.5 8 1.5 8])
        xlabel('SAUP P catch (MT)')
        ylabel('POEM P catch (MT)')
        title('Mean P catch')
        
        subplot(2,2,4)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,'--b'); hold on;
        plot(x,x2l,'--b'); hold on;
        plot(x,x5h,'--r'); hold on;
        plot(x,x5l,'--r'); hold on;
        plot(l10sD(keep),l10pD(keep),'.k','MarkerSize',25); hold on;
        text(2,7.5,['r = ' num2str(rD(k,b))])
        text(2,7.0,['RMSE = ' num2str(rmseD(k,b))])
        axis([1.5 8 1.5 8])
        xlabel('SAUP D catch (MT)')
        ylabel('POEM D catch (MT)')
        title('Mean D catch')
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished',harv,'_SAUP_log10catch_comp_types_LELC.png'])
        
        
        %% MAPS
        plme = NaN*ones(360,200);
        plmeF = NaN*ones(360,200);
        plmeP = NaN*ones(360,200);
        plmeD = NaN*ones(360,200);
        slme = NaN*ones(360,200);
        slmeF = NaN*ones(360,200);
        slmeP = NaN*ones(360,200);
        slmeD = NaN*ones(360,200);
        
        tlme = lme_mask_esm2m';
        [ni,nj]=size(geolon_t);
        
        for L=1:66
            lid = find(tlme==L);
            
            plme(lid) = plme_mcatch(L);
            plmeF(lid) = plme_Fmcatch(L);
            plmeP(lid) = plme_Pmcatch(L);
            plmeD(lid) = plme_Dmcatch(L);
            slme(lid) = slme_mcatch10(L);
            slmeF(lid) = Flme_mcatch10(L);
            slmeP(lid) = Plme_mcatch10(L);
            slmeD(lid) = Dlme_mcatch10(L);
        end
        
        % plot info
        geolon_t = double(geolon_t);
        geolat_t = double(geolat_t);
        plotminlat=-90; %Set these bounds for your data
        plotmaxlat=90;
        plotminlon=-280;
        plotmaxlon=80;
        latlim=[plotminlat plotmaxlat];
        lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
        % ENTER -100 TO MAP ORIGIN LONG
        
        land=-999*ones(ni,nj);
        land(grid(:,1))=NaN*ones(size(grid(:,1)));
        
        dlme = ((log10(plme)-log10(slme))./log10(slme)) .* 100;
        dlmeF = ((log10(plmeF)-log10(slmeF))./log10(slmeF)) .* 100;
        dlmeP = ((log10(plmeP)-log10(slmeP))./log10(slmeP)) .* 100;
        dlmeD = ((log10(plmeD)-log10(slmeD))./log10(slmeD)) .* 100;
        
        %% Catch
        
        % all
        figure(21)
        axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
            'Grid','off','FLineWidth',1,'origin',[0 -100 0])
        surfm(geolat_t,geolon_t,dlme)
        colormap(cmap_color_rb)
        load coast;                     %decent looking coastlines
        h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        caxis([-50 50]);
        hcb = colorbar('h');
        ylim(hcb,[-50 50])                   %Set color axis if needed
        set(gcf,'renderer','painters')
        title('Historic 1996-2005 difference from SAUP total annual catch (MT)')
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished_LME_SAUP_catch_diff.png'])
        
        %% F
        figure(22)
        axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
            'Grid','off','FLineWidth',1,'origin',[0 -100 0])
        surfm(geolat_t,geolon_t,dlmeF)
        colormap(cmap_color_rb)
        load coast;                     %decent looking coastlines
        h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        caxis([-50 50]);
        hcb = colorbar('h');
        ylim(hcb,[-50 50])                   %Set color axis if needed
        set(gcf,'renderer','painters')
        title('Historic 1996-2005 difference from SAUP total annual F catch (MT)')
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished_LME_SAUP_catch_diffF.png'])
        
        % all
        figure(23)
        axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
            'Grid','off','FLineWidth',1,'origin',[0 -100 0])
        surfm(geolat_t,geolon_t,dlmeP)
        colormap(cmap_color_rb)
        load coast;                     %decent looking coastlines
        h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        caxis([-50 50]);
        hcb = colorbar('h');
        ylim(hcb,[-50 50])                   %Set color axis if needed
        set(gcf,'renderer','painters')
        title('Historic 1996-2005 difference from SAUP total annual P catch (MT)')
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished_LME_SAUP_catch_diffP.png'])
        
        % all
        figure(24)
        axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
            'Grid','off','FLineWidth',1,'origin',[0 -100 0])
        surfm(geolat_t,geolon_t,dlmeD)
        colormap(cmap_color_rb)
        load coast;                     %decent looking coastlines
        h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        caxis([-50 50]);
        hcb = colorbar('h');
        ylim(hcb,[-50 50])                   %Set color axis if needed
        set(gcf,'renderer','painters')
        title('Historic 1996-2005 difference from SAUP total annual D catch (MT)')
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished_LME_SAUP_catch_diffD.png'])
        
    end
end

save([dp 'Hindcast_fished',harv,'_LME_SAUP_catch_comp_b_k.mat'],'rall','rF','rP',...
    'rD','rmse','rmseF','rmseP','rmseD')

%% Plots of r and RMSE
cfile2 = ['Hist_fished',harv,'_Dc_enc70_cmax-metab20_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100'];
        
kays2 = 0.0405:0.01:0.1016;
bees2 = 0.125:0.025:0.28;
[bgrid,kgrid]=meshgrid(bees2,kays2);
nc = length(kays);
ne = length(bees);

Rall  = NaN*ones(nc+1,ne+1);
RF  = NaN*ones(nc+1,ne+1);
RP  = NaN*ones(nc+1,ne+1);
RD  = NaN*ones(nc+1,ne+1);
Rall(1:nc,1:ne)=rall;
RF(1:nc,1:ne)=rF;
RP(1:nc,1:ne)=rP;
RD(1:nc,1:ne)=rD;

figure(1)
subplot(2,2,1)
pcolor(bgrid,kgrid,Rall)
colorbar
caxis([0.3 0.5])
set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays)
xlabel('Weight exponent')
ylabel('Temp exponent')
title('Corr all fish')

subplot(2,2,2)
pcolor(bgrid,kgrid,RF)
colorbar
caxis([-0.2 0.5])
set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays)
xlabel('Weight exponent')
ylabel('Temp exponent')
title('Corr F')

subplot(2,2,3)
pcolor(bgrid,kgrid,RP)
colorbar
caxis([0 0.7])
set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays)
xlabel('Weight exponent')
ylabel('Temp exponent')
title('Corr P')

subplot(2,2,4)
pcolor(bgrid,kgrid,RD)
colorbar
caxis([0 0.6])
set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays)
xlabel('Weight exponent')
ylabel('Temp exponent')
title('Corr D')
stamp(cfile2)
print('-dpng',[pp 'Hindcast_fished',harv,'_LME_SAUP_catch_corr_b_k'])


%%
Rmseall  = NaN*ones(nc+1,ne+1);
RmseF  = NaN*ones(nc+1,ne+1);
RmseP  = NaN*ones(nc+1,ne+1);
RmseD  = NaN*ones(nc+1,ne+1);
Rmseall(1:nc,1:ne)=rmse;
RmseF(1:nc,1:ne)=rmseF;
RmseP(1:nc,1:ne)=rmseP;
RmseD(1:nc,1:ne)=rmseD;

figure(2)
subplot(2,2,1)
pcolor(bgrid,kgrid,Rmseall)
colorbar
caxis([0.6 0.8])
set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays)
xlabel('Weight exponent')
ylabel('Temp exponent')
title('RMSE all fish')

subplot(2,2,2)
pcolor(bgrid,kgrid,RmseF)
colorbar
caxis([0.8 3])
set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays)
xlabel('Weight exponent')
ylabel('Temp exponent')
title('RMSE F')

subplot(2,2,3)
pcolor(bgrid,kgrid,RmseP)
colorbar
caxis([0.8 3])
set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays)
xlabel('Weight exponent')
ylabel('Temp exponent')
title('RMSE P')

subplot(2,2,4)
pcolor(bgrid,kgrid,RmseD)
colorbar
caxis([0.8 3])
set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays)
xlabel('Weight exponent')
ylabel('Temp exponent')
title('RMSE D')
stamp(cfile2)
print('-dpng',[pp 'Hindcast_fished',harv,'_LME_SAUP_catch_rmse_b_k'])

%%

corr_all = normalize(rall)+normalize(rF)+normalize(rP)+normalize(rD);
max(corr_all(:)) %k=0.09;b=0.2

rmse_all = normalize(rmse)+normalize(rmseF)+normalize(rmseP)+normalize(rmseD);
min(rmse_all(:)) %k=0.09;b=0.225


