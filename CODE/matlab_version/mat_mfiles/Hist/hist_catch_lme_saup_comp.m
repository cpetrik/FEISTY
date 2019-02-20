%POEM catch vs. SAUP catch by LME

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70_cmax-metab20_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00400';
harv = '03';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([dpath 'LME_hist_fished',harv,'_' cfile '.mat']);
load([spath 'LME_Catch_annual.mat']);
load([cpath 'hindcast_gridspec.mat'],'AREA_OCN','geolat_t','geolon_t');
load([cpath 'lme_mask_esm2m.mat']);
grid = csvread([cpath 'grid_csv.csv']);

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

%% 
%1955-2005 SAUP average
id = find(yr>1955 & yr<=2005);

slme_mcatch = nanmean(lme_catch(id,:));
slme_mcatch = slme_mcatch';


slme_mcatch10 = NaN*ones(size(slme_mcatch));
slme_mcatch20 = NaN*ones(size(slme_mcatch));
%Top 10 & 20 yrs SAUP
for i=1:66
    sort_lme_catch = sort(lme_catch(:,i),'descend');
    slme_mcatch10(i) = nanmean(sort_lme_catch(1:10));
    slme_mcatch20(i) = nanmean(sort_lme_catch(1:20));
end

%POEM LME biomass in MT
plme_mcatch = (mf_lme_mcatch+mp_lme_mcatch+md_lme_mcatch+...
    lp_lme_mcatch+ld_lme_mcatch) * 1e-6;

%Difference
diff_catch = plme_mcatch - slme_mcatch10;

code = [1:66]';
T = table(code,slme_mcatch10,plme_mcatch,'VariableNames',{'lme','saup','poem'});
writetable(T,[dpath 'LME_saup_catch_spinup_fished',harv,'.csv']);

%% Fit line
fraw = fit(slme_mcatch10,plme_mcatch,'poly1');
mraw = fraw.p1;
braw = fraw.p2;
yraw = mraw * slme_mcatch10 + braw;

l10s=log10(slme_mcatch10);
l10p=log10(plme_mcatch);
ii = find(l10s~=-Inf);
l10s = l10s(ii);
l10p = l10p(ii);

flog = fit(l10s,l10p,'poly1');
mlog = flog.p1;
blog = flog.p2;
ylog = mlog * l10s + blog;

%%
x=1:0.5:8;

figure(1)
plot(10.^x,10.^x,'--k'); hold on;
plot(slme_mcatch10,plme_mcatch,'.k','MarkerSize',25); hold on;
plot(slme_mcatch10,yraw,'r'); hold on;
axis([0 2e7 0 2e7])
xlabel('SAUP catch (MT) mean of top 10 years')
ylabel('POEM catch (MT) mean of Historic 1956-2005')
title('Mean catch')
print('-dpng',[ppath 'Historic_fished',harv,'_SAUP_catch_comp.png'])

%%
figure(2)
plot(x,x,'--k'); hold on;
plot(l10s,l10p,'.k','MarkerSize',25); hold on;
plot(l10s,ylog,'r'); hold on;
%axis([2 8 2 8])
xlabel('SAUP catch (MT) mean of top 10 years')
ylabel('POEM catch (MT) mean of Historic 1956-2005')
title('Mean catch')
print('-dpng',[ppath 'Historic_fished',harv,'_SAUP_log10catch_comp.png'])


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
plot(log10(slme_mcatch10(nwa)),log10(plme_mcatch(nwa)),'.c','MarkerSize',25); hold on;
axis([2 8 2 8])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('NW Atl 1956-2005 mean catch')
for i=1:length(nwa)
    text(log10(slme_mcatch10(nwa(i))),log10(plme_mcatch(nwa(i))),num2str(nwa(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(nea)),log10(plme_mcatch(nea)),'.c','MarkerSize',25); hold on;
axis([2 8 2 8])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('NE Atl')
for i=1:length(nea)
    text(log10(slme_mcatch10(nea(i))),log10(plme_mcatch(nea(i))),num2str(nea(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(nwp)),log10(plme_mcatch(nwp)),'.c','MarkerSize',25); hold on;
axis([2 8 2 8])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('NW Pac')
for i=1:length(nwp)
    text(log10(slme_mcatch10(nwp(i))),log10(plme_mcatch(nwp(i))),num2str(nwp(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(nep)),log10(plme_mcatch(nep)),'.c','MarkerSize',25); hold on;
axis([2 8 2 8])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('NE Pac')
for i=1:length(nep)
    text(log10(slme_mcatch10(nep(i))),log10(plme_mcatch(nep(i))),num2str(nep(i)),...
        'Color','k','HorizontalAlignment','center')
end
print('-dpng',[ppath 'Historic_fished',harv,'_SAUP10_log10catch_comp_AtlPac.png'])

%%
figure(6)
subplot(2,3,1)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(afr)),log10(plme_mcatch(afr)),'.c','MarkerSize',25); hold on;
axis([2 8 2 8])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('Africa 1956-2005 mean catch')
for i=1:length(afr)
    text(log10(slme_mcatch10(afr(i))),log10(plme_mcatch(afr(i))),num2str(afr(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,3,2)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(aus)),log10(plme_mcatch(aus)),'.c','MarkerSize',25); hold on;
axis([2 8 2 8])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('Australia')
for i=1:length(aus)
    text(log10(slme_mcatch10(aus(i))),log10(plme_mcatch(aus(i))),num2str(aus(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,3,3)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(ind)),log10(plme_mcatch(ind)),'.c','MarkerSize',25); hold on;
axis([2 8 2 8])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('India-Indonesia')
for i=1:length(ind)
    text(log10(slme_mcatch10(ind(i))),log10(plme_mcatch(ind(i))),num2str(ind(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,3,4)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(crb)),log10(plme_mcatch(crb)),'.c','MarkerSize',25); hold on;
axis([2 8 2 8])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('Caribbean')
for i=1:length(crb)
    text(log10(slme_mcatch10(crb(i))),log10(plme_mcatch(crb(i))),num2str(crb(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,3,5)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(sam)),log10(plme_mcatch(sam)),'.c','MarkerSize',25); hold on;
axis([2 8 2 8])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('S Amer')
for i=1:length(sam)
    text(log10(slme_mcatch10(sam(i))),log10(plme_mcatch(sam(i))),num2str(sam(i)),...
        'Color','k','HorizontalAlignment','center')
end
print('-dpng',[ppath 'Historic_fished',harv,'_SAUP10_log10catch_comp_Other.png'])

%% Drop Arctic, Antarctic, Hawaii, Australia

keep = [1:9,11:34,36:38,46:54,59:60,62,65];

figure(7)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(keep)),log10(plme_mcatch(keep)),'.k','MarkerSize',25); hold on;
axis([2 8 2 8])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('1956-2005 mean catch without Polar and Australia')
% for i=1:length(keep)
%     text(log10(slme_mcatch10(keep(i))),log10(plme_mcatch(keep(i))),num2str(keep(i)),...
%         'Color','k','HorizontalAlignment','center')
% end
print('-dpng',[ppath 'Historic_fished',harv,'_SAUP10_log10catch_comp_LELC.png'])

%% 
[rraw,praw] = corrcoef(slme_mcatch10,plme_mcatch)
[rlog,plog] = corrcoef(l10s,l10p)
[rlogO,plogO] = corrcoef(log10(slme_mcatch10(keep)),log10(plme_mcatch(keep)))

%%
figure(8)
plot(x,x,'--k'); hold on;
plot(l10s,l10p,'.k','MarkerSize',25); hold on;
axis([1 8 1 8])
xlabel('log10 SAUP catch (MT)')
ylabel('log10 POEM catch (MT)')
title('Mean catch')
text(2,6,['r = ' num2str(rlog(2,1))])
print('-dpng',[ppath 'Historic_fished',harv,'_SAUP_log10catch_comp_rp.png'])

%%
figure(9)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(keep)),log10(plme_mcatch(keep)),'.k','MarkerSize',25); hold on;
axis([1 8 1 8])
xlabel('log10 SAUP catch (MT)')
ylabel('log10 POEM catch (MT)')
text(2,6,['r = ' num2str(rlogO(2,1))])
title('Mean catch without Polar and Australia')
print('-dpng',[ppath 'Historic_fished',harv,'_SAUP_log10catch_comp_LELC_rp.png'])

%%
figure(10)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(keep)),log10(plme_mcatch(keep)),'.k','MarkerSize',25); hold on;
axis([4 8 4 8])
xlabel('log10 SAUP catch (MT)')
ylabel('log10 POEM catch (MT)')
text(5,7.5,['r = ' num2str(rlogO(2,1))])
title('Mean catch without Australia and polar regions')
print('-dpng',[ppath 'Historic_fished',harv,'_SAUP_log10catch_comp_LELC_rp_crop.png'])

%%
figure(11)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10),log10(plme_mcatch),'.k','MarkerSize',25); hold on;
%axis([1 8 1 8])
xlabel('log10 SAUP catch (MT)')
ylabel('log10 POEM catch (MT)')
title('Mean catch')
text(2,6,['r = ' num2str(rraw(2,1))])
print('-dpng',[ppath 'Historic_fished',harv,'_SAUP_log10catch_raw_rp.png'])

%% MAPS
plme = NaN*ones(360,200);
slme = NaN*ones(360,200);

tlme = lme_mask_esm2m';
[ni,nj]=size(geolon_t);

for L=1:66
    lid = find(tlme==L);
    
    plme(lid) = plme_mcatch(L);
    slme(lid) = slme_mcatch(L);
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

dlme = ((plme-slme)./slme) .* 100;

%% Catch

% all
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,dlme)
colormap(cmap_color_rb)            
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
hcb = colorbar('h');
ylim(hcb,[-100 100])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic 1956-2005 difference from SAUP total annual catch (MT)')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_LME_SAUP_catch_diff.png'])


