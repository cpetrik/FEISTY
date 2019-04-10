function clim_catch_lme_saup_Fcorr_stock_LMEid(lme_Fmcatch,ppath,simname)

%FEISTY catch vs. SAUP catch by LME
%Use same methods as Stock et al. 2017 to reduce SAUP dataset

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


%% Assign a color to each LME based on temp
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
load([cpath 'LME_clim_temp_zoop_det.mat'],'lme_ptemp');

tmap=cmocean('thermal',66);
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

%% SAUP
% MT/km2
l10sF=log10(Flme_mcatch10+eps);

%% FEISTY LME biomass in MT
plme_Fmcatch = lme_Fmcatch * 1e-6;
% MT/km2
plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;

l10pF=log10(plme_Fmcatch);

%% F Major countries
figure(1)
subplot(3,3,1)
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
%scatter(l10sF(keep),l10pF(keep),20,lme_ptemp(keep,1),'filled'); hold on;
%cmocean('thermal');
axis([-6 2 -6 2])
title('F NW Atl mean catch')

subplot(3,3,2)
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
title('F NE Atl')

subplot(3,3,3)
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
title('F NW Pac')

subplot(3,3,4)
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
ylabel('log10 FEISTY catch (MT km^-^2)')
title('F NE Pac')

% F Minor countries
subplot(3,3,5)
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
axis([-6 2 -6 2])
title('F Africa mean catch')

subplot(3,3,6)
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
axis([-6 2 -6 2])
title('F Australia')

subplot(3,3,7)
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
axis([-6 2 -6 2])
title('F India-Indonesia')

subplot(3,3,8)
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
axis([-6 2 -6 2])
xlabel('log10 SAUP catch (MT km^-^2) top 10 yrs')
title('F Caribbean')

subplot(3,3,9)
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
axis([-6 2 -6 2])
title('F S Amer')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_SAUP_comp_Stock_compF_',simname,'.png'])


end
