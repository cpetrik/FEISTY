%POEM catch vs. SAUP catch by LME

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Diff_Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE10_CC100_RE00100';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([dpath 'LME_spinup_fished_' cfile '.mat']);
load([cpath 'hindcast_gridspec.mat'],'AREA_OCN','geolat_t','geolon_t');
load([cpath 'LME_clim_temp.mat']);
load([cpath 'lme_mask_esm2m.mat']);
grid = csvread([cpath 'grid_csv.csv']);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);
lme_area_km2 = lme_area * 1e-6; 

%% SAUP catch in MT
%use weighted catches
load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
    'Flme_wcatch','Dlme_wcatch','Plme_wcatch');
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

%% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);

%non-LELC
keep = [1:9,11:34,36:38,46:54,59:60,62,65];

%Lines
x=-6:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

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
ylabel('POEM mean of Spinup')
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
print('-dpng',[ppath 'Spinup_fished_SAUP_log10catch_comp_types_temp.png'])



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


figure(2)
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
xlabel('SAUP mean of top 10 years')
ylabel('POEM total mean of Spinup')
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
print('-dpng',[ppath 'Spinup_fished_SAUP_log10catch_comp_types_LELC_temp.png'])


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

figure(5)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(l10s(nwa),l10p(nwa),'.c','MarkerSize',25); hold on;
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('NW Atl 1951-2000 mean catch')
for i=1:length(nwa)
    text(l10s(nwa(i)),l10p(nwa(i)),num2str(nwa(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(l10s(nea),l10p(nea),'.c','MarkerSize',25); hold on;
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('NE Atl')
for i=1:length(nea)
    text(l10s(nea(i)),l10p(nea(i)),num2str(nea(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(l10s(nwp),l10p(nwp),'.c','MarkerSize',25); hold on;
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('NW Pac')
for i=1:length(nwp)
    text(l10s(nwp(i)),l10p(nwp(i)),num2str(nwp(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(l10s(nep),l10p(nep),'.c','MarkerSize',25); hold on;
axis([-4 2 -4 2])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('NE Pac')
for i=1:length(nep)
    text(l10s(nep(i)),l10p(nep(i)),num2str(nep(i)),...
        'Color','k','HorizontalAlignment','center')
end
print('-dpng',[ppath 'Spinup_fished_SAUP10_log10catch_comp_AtlPac.png'])

%%
figure(6)
subplot(2,3,1)
plot(x,x,'--k'); hold on;
plot(l10s(afr),l10p(afr),'.c','MarkerSize',25); hold on;
axis([-2 2 -2 2])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('Africa 1951-2000 mean catch')
for i=1:length(afr)
    text(l10s(afr(i)),l10p(afr(i)),num2str(afr(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,3,2)
plot(x,x,'--k'); hold on;
plot(l10s(aus),l10p(aus),'.c','MarkerSize',25); hold on;
axis([-2 2 -2 2])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('Australia')
for i=1:length(aus)
    text(l10s(aus(i)),l10p(aus(i)),num2str(aus(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,3,3)
plot(x,x,'--k'); hold on;
plot(l10s(ind),l10p(ind),'.c','MarkerSize',25); hold on;
axis([-2 2 -2 2])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('India-Indonesia')
for i=1:length(ind)
    text(l10s(ind(i)),l10p(ind(i)),num2str(ind(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,3,4)
plot(x,x,'--k'); hold on;
plot(l10s(crb),l10p(crb),'.c','MarkerSize',25); hold on;
axis([-2 2 -2 2])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('Caribbean')
for i=1:length(crb)
    text(l10s(crb(i)),l10p(crb(i)),num2str(crb(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,3,5)
plot(x,x,'--k'); hold on;
plot(l10s(sam),l10p(sam),'.c','MarkerSize',25); hold on;
axis([-2 2 -2 2])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('S Amer')
for i=1:length(sam)
    text(l10s(sam(i)),l10p(sam(i)),num2str(sam(i)),...
        'Color','k','HorizontalAlignment','center')
end
print('-dpng',[ppath 'Spinup_fished_SAUP10_log10catch_comp_Other.png'])

