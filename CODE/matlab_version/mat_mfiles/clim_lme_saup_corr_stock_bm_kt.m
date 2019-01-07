%FEISTY catch vs. SAUP catch by LME
%Use same methods as Stock et al. 2017 to reduce SAUP dataset

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'LME_clim_temp_zoop_det.mat']);

%use weighted catches
load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
    'Flme_wcatch','Dlme_wcatch','Plme_wcatch');
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')

AREA_OCN = max(area,1);

% FEISTY file info
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = [pp cfile '/param_sens/'];

ktemp = 0.0705:0.01:0.1105;
bees = -0.15:-0.025:-0.25;

load(['/Volumes/GFDL/NC/Matlab_new_size/',cfile,...
    '/param_sens/LME_Climatol_All_fish03_bm_kt_search.mat']);

lme_area_km2 = lme_area * 1e-6;
lme_area_km2_loop = NaN*ones(length(bees),length(ktemp),66);
for k=1:length(ktemp)
    for j=1:length(bees)
        lme_area_km2_loop(j,k,:) = lme_area_km2;
    end
end

%% DvD on grid
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/DanielVD_PelDem/Colleen_modeledfish_LME.mat')
dlme_Pfrac = NaN*ones(180,360);
for L=1:63
    lid = find(tlme==L);
    dlme_Pfrac(lid) = FracLP(L);
end

did=[1:61,63];

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

sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

%% plot info
x=-8:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

cmYOR=cbrewer('seq','YlOrRd',28);

%% FEISTY LME biomass in MT
plme_mcatch = (lme_mcatch_mf+lme_mcatch_mp+lme_mcatch_lp+lme_mcatch_md+...
    lme_mcatch_ld) * 1e-6;
plme_Fmcatch = (lme_mcatch_mf) * 1e-6;
plme_Pmcatch = (lme_mcatch_mp+lme_mcatch_lp) * 1e-6;
plme_Dmcatch = (lme_mcatch_md+lme_mcatch_ld) * 1e-6;
plme_Bmbio = lme_mbio_b * 1e-6;
% MT/km2
plme_mcatch = plme_mcatch ./ lme_area_km2_loop;
plme_Fmcatch = plme_Fmcatch ./ lme_area_km2_loop;
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2_loop;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2_loop;
plme_Bmbio = plme_Bmbio ./ lme_area_km2_loop;

pFracPD = plme_Pmcatch ./ (plme_Pmcatch + plme_Dmcatch);

l10p=log10(plme_mcatch);
l10pF=log10(plme_Fmcatch);
l10pP=log10(plme_Pmcatch);
l10pD=log10(plme_Dmcatch);

%% Drop Arctic, Antarctic, Hawaii, Australia -------------------------
rall = NaN*ones(length(bees),length(ktemp));
rF = rall;
rP = rall;
rD = rall;
rPD = rall;
rVD = rall;
pall = rall;
pF = rall;
pP = rall;
pD = rall;
pPD = rall;
pVD = rall;
rmse = rall;
rmseF = rall;
rmseP = rall;
rmseD = rall;
rmsePD = rall;
rmseVD = rall;
for k=1:length(ktemp)
    for j=1:length(bees)
        % Stats
        %r
        [rall(j,k),pall(j,k)]=corr(l10s(keep),squeeze(l10p(j,k,keep)));
        [rF(j,k),pF(j,k)]=corr(l10sF(keep),squeeze(l10pF(j,k,keep)));
        [rP(j,k),pP(j,k)]=corr(l10sP(keep),squeeze(l10pP(j,k,keep)));
        [rD(j,k),pD(j,k)]=corr(l10sD(keep),squeeze(l10pD(j,k,keep)));
        [rPD(j,k),pPD(j,k)]=corr(sFracPD(keep),squeeze(pFracPD(j,k,keep)));
        [rVD(j,k),pVD(j,k)]=corr(FracLP(did),squeeze(pFracPD(j,k,did)));
        
        %root mean square error
        o=l10s(keep);
        p=squeeze(l10p(j,k,keep));
        n = length(o);
        num=nansum((p-o).^2);
        rmse(j,k) = sqrt(num/n);
        
        o=l10sF(keep);
        p=squeeze(l10pF(j,k,keep));
        n = length(o);
        num=nansum((p-o).^2);
        rmseF(j,k) = sqrt(num/n);
        
        o=l10sP(keep);
        p=squeeze(l10pP(j,k,keep));
        n = length(o);
        num=nansum((p-o).^2);
        rmseP(j,k) = sqrt(num/n);
        
        o=l10sD(keep);
        p=squeeze(l10pD(j,k,keep));
        n = length(o);
        num=nansum((p-o).^2);
        rmseD(j,k) = sqrt(num/n);
        
        o=sFracPD(keep);
        p=squeeze(pFracPD(j,k,keep));
        n = length(o);
        num=nansum((p-o).^2);
        rmsePD(j,k) = sqrt(num/n);
        
        o=FracLP(did);
        p=squeeze(pFracPD(j,k,did));
        n = length(o);
        num=nansum((p-o).^2);
        rmseVD(j,k) = sqrt(num/n);
        
    end
end

%% Plots
n=0;
for j=1:5
    for k=1:length(ktemp)
        n=n+1;
        % For ms
        f1=figure(1);
        subplot(5,5,n)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,':b'); hold on;
        plot(x,x2l,':b'); hold on;
        plot(x,x5h,':r'); hold on;
        plot(x,x5l,':r'); hold on;
        scatter(l10sF(keep),squeeze(l10pF(j,k,keep)),20,lme_ptemp(keep,1),'filled'); hold on;
        cmocean('thermal');
        colorbar('Position',[0.92 0.4 0.025 0.3],'orientation','vertical')
        text(-5.5,1.5,['r = ' sprintf('%2.2f',rF(j,k))])
        text(-5.5,-5.5,['E = ' sprintf('%2.2f',rmseF(j,k))])
        axis([-6 2 -6 2])
        if(k==1)
            ylabel(num2str(bees(j)),'FontWeight','bold')
        end
        if (n==11)
            ylabel({'FEISTY';num2str(bees(j))},'FontWeight','bold')
        end
        if (n==23)
            xlabel('SAU','FontWeight','bold')
        end
        if(j==5)
            title(num2str(ktemp(k)),'FontWeight','bold')
        end
        if (n==3)
            title({'Forage Fishes';num2str(ktemp(k))})
            stamp(cfile)
        end
        
        f2=figure(2);
        subplot(5,5,n)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,':b'); hold on;
        plot(x,x2l,':b'); hold on;
        plot(x,x5h,':r'); hold on;
        plot(x,x5l,':r'); hold on;
        scatter(l10sP(keep),squeeze(l10pP(j,k,keep)),20,lme_ptemp(keep,1),'filled'); hold on;
        cmocean('thermal');
        colorbar('Position',[0.92 0.4 0.025 0.3],'orientation','vertical')
        text(-5.5,1.5,['r = ' sprintf('%2.2f',rP(j,k))])
        text(-5.5,-5.5,['E = ' sprintf('%2.2f',rmseP(j,k))])
        axis([-6 2 -6 2])
        if(k==1)
            ylabel(num2str(bees(j)),'FontWeight','bold')
        end
        if (n==11)
            ylabel({'FEISTY';num2str(bees(j))},'FontWeight','bold')
        end
        if (n==23)
            xlabel('SAU','FontWeight','bold')
        end
        if(j==5)
            title(num2str(ktemp(k)),'FontWeight','bold')
        end
        if (n==3)
            title({'Large Pelagics';num2str(ktemp(k))})
            stamp(cfile)
        end
        
        f3=figure(3);
        subplot(5,5,n)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,':b'); hold on;
        plot(x,x2l,':b'); hold on;
        plot(x,x5h,':r'); hold on;
        plot(x,x5l,':r'); hold on;
        scatter(l10sD(keep),squeeze(l10pD(j,k,keep)),20,lme_ptemp(keep,1),'filled'); hold on;
        cmocean('thermal');
        colorbar('Position',[0.92 0.4 0.025 0.3],'orientation','vertical')
        text(-1.75,0.75,['r = ' sprintf('%2.2f',rD(j,k))])
        text(-1.75,-1.75,['E = ' sprintf('%2.2f',rmseD(j,k))])
        axis([-2 1 -2 1])
        if(k==1)
            ylabel(num2str(bees(j)),'FontWeight','bold')
        end
        if (n==11)
            ylabel({'FEISTY';num2str(bees(j))},'FontWeight','bold')
        end
        if (n==23)
            xlabel('SAU','FontWeight','bold')
        end
        if(j==5)
            title(num2str(ktemp(k)),'FontWeight','bold')
        end
        if (n==3)
            title({'Demersals';num2str(ktemp(k))})
            stamp(cfile)
        end
        
        f4=figure(4);
        subplot(5,5,n)
        plot(x,x,'--k'); hold on;
        plot(x,x2h,':b'); hold on;
        plot(x,x2l,':b'); hold on;
        plot(x,x5h,':r'); hold on;
        plot(x,x5l,':r'); hold on;
        scatter(l10s(keep),squeeze(l10p(j,k,keep)),20,lme_ptemp(keep,1),'filled'); hold on;
        cmocean('thermal');
        colorbar('Position',[0.92 0.4 0.025 0.3],'orientation','vertical')
        text(-1.75,1.75,['r = ' sprintf('%2.2f',rall(j,k))])
        text(-1.75,-1.75,['E = ' sprintf('%2.2f',rmse(j,k))])
        axis([-2 2 -2 2])
        if(k==1)
            ylabel(num2str(bees(j)),'FontWeight','bold')
        end
        if (n==11)
            ylabel({'FEISTY';num2str(bees(j))},'FontWeight','bold')
        end
        if (n==23)
            xlabel('SAU','FontWeight','bold')
        end
        if(j==5)
            title(num2str(ktemp(k)),'FontWeight','bold')
        end
        if (n==3)
            title({'All Fishes';num2str(ktemp(k))})
            stamp(cfile)
        end
        
        %SAUP frac
        f5 = figure(5);
        subplot(5,5,n)
        plot(x,x,'--k'); hold on;
        scatter(sFracPD(keep),squeeze(pFracPD(j,k,keep)),20,lme_ptemp(keep,1),'filled'); hold on;
        cmocean('thermal');
        colorbar('Position',[0.92 0.4 0.025 0.3],'orientation','vertical')
        text(0.1,0.95,['r = ' sprintf('%2.2f',rPD(j,k))])
        text(0.5,0.1,['E = ' sprintf('%2.2f',rmsePD(j,k))])
        axis([0 1.05 0 1.05])
        if(k==1)
            ylabel(num2str(bees(j)),'FontWeight','bold')
        end
        if (n==11)
            ylabel({'FEISTY';num2str(bees(j))},'FontWeight','bold')
        end
        if (n==23)
            xlabel('SAU','FontWeight','bold')
        end
        if(j==5)
            title(num2str(ktemp(k)),'FontWeight','bold')
        end
        if (n==3)
            title({'P/(P+D)';num2str(ktemp(k))})
            stamp(cfile)
        end
        
        %DvD Corr
        f6 = figure(6);
        subplot(5,5,n)
        plot(x,x,'--k'); hold on;
        scatter(FracLP(did),squeeze(pFracPD(j,k,did)),20,lme_ptemp(did,1),'filled'); hold on;
        cmocean('thermal');
        colorbar('Position',[0.92 0.4 0.025 0.3],'orientation','vertical')
        text(0.1,0.95,['r = ' sprintf('%2.2f',rVD(j,k))])
        text(0.5,0.1,['E = ' sprintf('%2.2f',rmseVD(j,k))])
        axis([0 1.05 0 1.05])
        if(k==1)
            ylabel(num2str(bees(j)),'FontWeight','bold')
        end
        if (n==11)
            ylabel({'FEISTY';num2str(bees(j))},'FontWeight','bold')
        end
        if (n==23)
            xlabel('vanD','FontWeight','bold')
        end
        if(j==5)
            title(num2str(ktemp(k)),'FontWeight','bold')
        end
        if (n==3)
            title({'P/(P+D)';num2str(ktemp(k))})
            stamp(cfile)
        end
        
    end
end
%%
print(f1,'-dpng',[ppath 'Clim_fish03_SAUP_comp_Stock_LELC_bm_kt_F.png'])
print(f2,'-dpng',[ppath 'Clim_fish03_SAUP_comp_Stock_LELC_bm_kt_P.png'])
print(f3,'-dpng',[ppath 'Clim_fish03_SAUP_comp_Stock_LELC_bm_kt_D.png'])
print(f4,'-dpng',[ppath 'Clim_fish03_SAUP_comp_Stock_LELC_bm_kt_All.png'])
print(f5,'-dpng',[ppath 'Clim_fish03_SAUP_comp_Stock_LELC_bm_kt_PD.png'])
print(f6,'-dpng',[ppath 'Clim_fish03_SAUP_comp_Stock_LELC_bm_kt_PD_DvD.png'])

%% pcolor
%ktemp = 0.0705:0.01:0.1105;
%bees = 0.15:0.025:0.25;
rbmp = -0.25:0.025:-0.15;
jays = [bees -0.275];
ays = [ktemp 0.1205];
[agrid,jgrid]=meshgrid(ays,jays);

nj = length(bees);
nk = length(ktemp);

rF2 = NaN*ones(nj+1,nk+1);
rP2 = NaN*ones(nj+1,nk+1);
rD2 = NaN*ones(nj+1,nk+1);
r2 = NaN*ones(nj+1,nk+1);
rPD2 = NaN*ones(nj+1,nk+1);
rVD2 = NaN*ones(nj+1,nk+1);
eF2 = NaN*ones(nj+1,nk+1);
eP2 = NaN*ones(nj+1,nk+1);
eD2 = NaN*ones(nj+1,nk+1);
e2 = NaN*ones(nj+1,nk+1);
ePD2 = NaN*ones(nj+1,nk+1);
eVD2 = NaN*ones(nj+1,nk+1);

rF2(1:nj,1:nk,:)=rF;
rP2(1:nj,1:nk,:)=rP;
rD2(1:nj,1:nk,:)=rD;
r2(1:nj,1:nk,:)=rall;
rPD2(1:nj,1:nk,:)=rPD;
rVD2(1:nj,1:nk,:)=rVD;
eF2(1:nj,1:nk,:)=rmseF;
eP2(1:nj,1:nk,:)=rmseP;
eD2(1:nj,1:nk,:)=rmseD;
e2(1:nj,1:nk,:)=rmse;
ePD2(1:nj,1:nk,:)=rmsePD;
eVD2(1:nj,1:nk,:)=rmseVD;

%%
figure(7);
subplot(3,3,1)
pcolor(agrid,jgrid,rF2)
colormap(cmYOR)
caxis([0 1])
ylabel('b_M')
title('F')
set(gca,'YTick',rbmp(1:2:end)+0.0125,'YTickLabel',rbmp(1:2:end),...
        'XTick',ktemp(1:2:end)+0.005,'XTickLabel',[0.07:0.02:0.11])

subplot(3,3,2)
pcolor(agrid,jgrid,rP2)
colormap(cmYOR)
caxis([0 1])
title('P')
set(gca,'YTick',bees(1:2:end)+0.0125,'YTickLabel',bees(1:2:end),...
        'XTick',ktemp(1:2:end)+0.005,'XTickLabel',[0.07:0.02:0.11])

subplot(3,3,3)
pcolor(agrid,jgrid,rD2)
colormap(cmYOR)
caxis([0 1])
title('D')
set(gca,'YTick',bees(1:2:end)+0.0125,'YTickLabel',bees(1:2:end),...
        'XTick',ktemp(1:2:end)+0.005,'XTickLabel',[0.07:0.02:0.11])

subplot(3,3,4)
pcolor(agrid,jgrid,r2)
colormap(cmYOR)
caxis([0 1])
xlabel('k_M')
ylabel('b_M')
title('All')
set(gca,'YTick',bees(1:2:end)+0.0125,'YTickLabel',bees(1:2:end),...
        'XTick',ktemp(1:2:end)+0.005,'XTickLabel',[0.07:0.02:0.11])

subplot(3,3,5)
pcolor(agrid,jgrid,rPD2)
colormap(cmYOR)
caxis([0 1])
xlabel('k_M')
title('SAU P/(P+D)')
set(gca,'YTick',bees(1:2:end)+0.0125,'YTickLabel',bees(1:2:end),...
        'XTick',ktemp(1:2:end)+0.005,'XTickLabel',[0.07:0.02:0.11])

subplot(3,3,6)
pcolor(agrid,jgrid,rVD2)
colorbar('Position',[0.375 0.325 0.3 0.03],'orientation','horizontal')
colormap(cmYOR)
caxis([0 1])
xlabel('k_M')
title('vanD P/(P+D)')
set(gca,'YTick',bees(1:2:end)+0.0125,'YTickLabel',bees(1:2:end),...
        'XTick',ktemp(1:2:end)+0.005,'XTickLabel',[0.07:0.02:0.11])
print('-dpng',[ppath 'Clim_fish03_SAUP_comp_Stock_LELC_bm_kt_pcolor_corrs.png'])    


figure(8);
subplot(3,3,1)
pcolor(agrid,jgrid,eF2)
colormap(cmYOR)
caxis([0 1])
ylabel('b_M')
title('F')
set(gca,'YTick',bees(1:2:end)+0.0125,'YTickLabel',bees(1:2:end),...
        'XTick',ktemp(1:2:end)+0.005,'XTickLabel',[0.07:0.02:0.11])

subplot(3,3,2)
pcolor(agrid,jgrid,eP2)
colormap(cmYOR)
caxis([0 1])
title('P')
set(gca,'YTick',bees(1:2:end)+0.0125,'YTickLabel',bees(1:2:end),...
        'XTick',ktemp(1:2:end)+0.005,'XTickLabel',[0.07:0.02:0.11])

subplot(3,3,3)
pcolor(agrid,jgrid,eD2)
colormap(cmYOR)
caxis([0 1])
title('D')
set(gca,'YTick',bees(1:2:end)+0.0125,'YTickLabel',bees(1:2:end),...
        'XTick',ktemp(1:2:end)+0.005,'XTickLabel',[0.07:0.02:0.11])

subplot(3,3,4)
pcolor(agrid,jgrid,e2)
colormap(cmYOR)
caxis([0 1])
xlabel('k_M')
ylabel('b_M')
title('All')
set(gca,'YTick',bees(1:2:end)+0.0125,'YTickLabel',bees(1:2:end),...
        'XTick',ktemp(1:2:end)+0.005,'XTickLabel',[0.07:0.02:0.11])

subplot(3,3,5)
pcolor(agrid,jgrid,ePD2)
colormap(cmYOR)
caxis([0 1])
xlabel('k_M')
title('SAU P/(P+D)')
set(gca,'YTick',bees(1:2:end)+0.0125,'YTickLabel',bees(1:2:end),...
        'XTick',ktemp(1:2:end)+0.005,'XTickLabel',[0.07:0.02:0.11])

subplot(3,3,6)
pcolor(agrid,jgrid,eVD2)
colorbar('Position',[0.375 0.325 0.3 0.03],'orientation','horizontal')
colormap(cmYOR)
caxis([0 1])
xlabel('k_M')
title('vanD P/(P+D)')
set(gca,'YTick',bees(1:2:end)+0.0125,'YTickLabel',bees(1:2:end),...
        'XTick',ktemp(1:2:end)+0.005,'XTickLabel',[0.07:0.02:0.11])
print('-dpng',[ppath 'Clim_fish03_SAUP_comp_Stock_LELC_bm_kt_pcolor_rmse.png'])    


