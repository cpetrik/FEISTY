%POEM catch vs. SAUP catch by LME
%Use same methods as Stock et al. 2017 to reduce SAUP dataset

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
% load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
% load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
% load([cpath 'esm26_area_1deg.mat']);
load([cpath 'LME_clim_temp_zoop_det.mat']);

%use weighted catches
load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
    'Flme_wcatch','Dlme_wcatch','Plme_wcatch');

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')

% "lfile" never changes, has lme areas
lfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
lpath = ['/Volumes/GFDL/NC/Matlab_new_size/' lfile '/'];
load([lpath 'LME_clim_fished03_' lfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

% POEM file info
frate = 0.3;
tfish = num2str(100+int64(10*frate));

cfile = 'Dc_enc70-b200_cm25_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat'],'lme_mcatch','lme_mbio','lme_sbio');
%load([dpath 'LME_clim_',harv,'_loop_' cfile '.mat'],'lme_mcatch');


%% Assign a color to each LME based on temp
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
nwa = [6:9,18];
nea = [19:22,24:26,59:60];
nep = [1:4,11,54,65];
nwp = [47,49:53,56];
sam = 13:17;
afr = [27:31];
aus = 46;
ind = [32,34,36:38];
crb = [5,12];

x=-6:0.5:8;
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

%% POEM LME biomass in MT
plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
plme_Bmbio = lme_mbio(:,9) * 1e-6;
plme_Bsbio = lme_sbio(:,9) * 1e-6;
% MT/km2
plme_mcatch = plme_mcatch ./ lme_area_km2;
plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;
plme_Bmbio = plme_Bmbio ./ lme_area_km2;
plme_Bsbio = plme_Bsbio ./ lme_area_km2;

l10p=log10(plme_mcatch);
l10pF=log10(plme_Fmcatch);
l10pP=log10(plme_Pmcatch);
l10pD=log10(plme_Dmcatch);

% Plots -------------------------------------------
%% F Major countries
figure(1)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(nwa)
    lme=nwa(i);
    plot(l10sF(lme),l10pF(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sF(lme),l10pF(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('F NW Atl mean catch')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(nea)
    lme=nea(i);
    plot(l10sF(lme),l10pF(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sF(lme),l10pF(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('F NE Atl')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(nwp)
    lme=nwp(i);
    plot(l10sF(lme),l10pF(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sF(lme),l10pF(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('F NW Pac')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(nep)
    lme=nep(i);
    plot(l10sF(lme),l10pF(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sF(lme),l10pF(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('F NE Pac')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_fished',harv,'_SAUP10_log10catch_compF_AtlPac.png'])

%% F Minor countries
figure(2)
subplot(2,3,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(afr)
    lme=afr(i);
    plot(l10sF(lme),l10pF(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sF(lme),l10pF(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('F Africa mean catch')

subplot(2,3,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(aus)
    lme=aus(i);
    plot(l10sF(lme),l10pF(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sF(lme),l10pF(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('F Australia')

subplot(2,3,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(ind)
    lme=ind(i);
    plot(l10sF(lme),l10pF(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sF(lme),l10pF(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('F India-Indonesia')

subplot(2,3,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(crb)
    lme=crb(i);
    plot(l10sF(lme),l10pF(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sF(lme),l10pF(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('F Caribbean')

subplot(2,3,5)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(sam)
    lme=sam(i);
    plot(l10sF(lme),l10pF(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sF(lme),l10pF(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('F S Amer')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_fished',harv,'_SAUP10_log10catch_compF_Other.png'])


%% P major countries
figure(3)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(nwa)
    lme=nwa(i);
    plot(l10sP(lme),l10pP(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sP(lme),l10pP(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('P NW Atl mean catch')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(nea)
    lme=nea(i);
    plot(l10sP(lme),l10pP(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sP(lme),l10pP(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('P NE Atl')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(nwp)
    lme=nwp(i);
    plot(l10sP(lme),l10pP(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sP(lme),l10pP(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('P NW Pac')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(nep)
    lme=nep(i);
    plot(l10sP(lme),l10pP(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sP(lme),l10pP(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('P NE Pac')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_fished',harv,'_SAUP10_log10catch_compP_AtlPac.png'])

%% P Minor countries
figure(4)
subplot(2,3,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(afr)
    lme=afr(i);
    plot(l10sP(lme),l10pP(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sP(lme),l10pP(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('P Africa mean catch')

subplot(2,3,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(aus)
    lme=aus(i);
    plot(l10sP(lme),l10pP(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sP(lme),l10pP(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('P Australia')

subplot(2,3,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(ind)
    lme=ind(i);
    plot(l10sP(lme),l10pP(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sP(lme),l10pP(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('P India-Indonesia')

subplot(2,3,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(crb)
    lme=crb(i);
    plot(l10sP(lme),l10pP(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sP(lme),l10pP(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('P Caribbean')

subplot(2,3,5)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(sam)
    lme=sam(i);
    plot(l10sP(lme),l10pP(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sP(lme),l10pP(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('P S Amer')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_fished',harv,'_SAUP10_log10catch_compP_Other.png'])

%% D major countries
figure(5)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(nwa)
    lme=nwa(i);
    plot(l10sD(lme),l10pD(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sD(lme),l10pD(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('D NW Atl mean catch')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(nea)
    lme=nea(i);
    plot(l10sD(lme),l10pD(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sD(lme),l10pD(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('D NE Atl')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(nwp)
    lme=nwp(i);
    plot(l10sD(lme),l10pD(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sD(lme),l10pD(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('D NW Pac')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(nep)
    lme=nep(i);
    plot(l10sD(lme),l10pD(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sD(lme),l10pD(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('D NE Pac')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_fished',harv,'_SAUP10_log10catch_compD_AtlPac.png'])

%% D Minor countries
figure(6)
subplot(2,3,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(afr)
    lme=afr(i);
    plot(l10sD(lme),l10pD(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sD(lme),l10pD(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('D Africa mean catch')

subplot(2,3,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(aus)
    lme=aus(i);
    plot(l10sD(lme),l10pD(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sD(lme),l10pD(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('D Australia')

subplot(2,3,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(ind)
    lme=ind(i);
    plot(l10sD(lme),l10pD(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sD(lme),l10pD(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('D India-Indonesia')

subplot(2,3,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(crb)
    lme=crb(i);
    plot(l10sD(lme),l10pD(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sD(lme),l10pD(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('D Caribbean')

subplot(2,3,5)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(sam)
    lme=sam(i);
    plot(l10sD(lme),l10pD(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(l10sD(lme),l10pD(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center')
end
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('D S Amer')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_fished',harv,'_SAUP10_log10catch_compD_Other.png'])
