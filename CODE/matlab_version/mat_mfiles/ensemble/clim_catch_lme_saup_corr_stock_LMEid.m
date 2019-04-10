%POEM catch vs. SAUP catch by LME
%Use same methods as Stock et al. 2017 to reduce SAUP dataset

clear all
close all

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/SAUP_Stock_top10.mat');

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
dpath = ['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

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
% MT/km2
sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

%% FEISTY  LME biomass in MT
plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
% MT/km2
plme_mcatch = plme_mcatch ./ lme_area_km2;
plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;

pFracPD = plme_Pmcatch ./ (plme_Pmcatch + plme_Dmcatch);

l10p=log10(plme_mcatch+eps);
l10pF=log10(plme_Fmcatch+eps);
l10pP=log10(plme_Pmcatch+eps);
l10pD=log10(plme_Dmcatch+eps);

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
axis([-6 2 -6 2])
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
axis([-6 2 -6 2])
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
axis([-6 2 -6 2])
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
axis([-6 2 -6 2])
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
axis([-3 2 -3 2])
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
axis([-3 2 -3 2])
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
axis([-3 2 -3 2])
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
axis([-3 2 -3 2])
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
axis([-3 2 -3 2])
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
axis([-6 1 -6 1])
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
axis([-6 1 -6 1])
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
axis([-6 1 -6 1])
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
axis([-6 1 -6 1])
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
axis([-5 2 -5 2])
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
axis([-2 1 -2 1])
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
axis([-3 1 -3 1])
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
axis([-3 0 -3 0])
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
axis([-3 2 -3 2])
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
axis([-2 2 -2 2])
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
axis([-2 2 -2 2])
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
axis([-2 2 -2 2])
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
axis([-2 2 -2 2])
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
axis([-2 1 -2 1])
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
axis([-2 1 -2 1])
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
axis([-2 1 -2 1])
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
axis([-2 1 -2 1])
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
axis([-2 1 -2 1])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
ylabel('log10 POEM catch (MT km^-^2)')
title('D S Amer')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_fished',harv,'_SAUP10_log10catch_compD_Other.png'])
