% Plot effective TEs at LME scale
% Climatology
% 150 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);

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
% ENTER -100 TO MAP ORIGIN LONG

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));


%%
AREA_OCN = max(area,1);

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'TEeff_Climatol_All_fish03_' cfile '.mat']);

%% Calc LMEs
tlme = lme_mask_onedeg;

lme_te = NaN*ones(66,2);
for L=1:66
    lid = find(tlme==L);
    %TEeff
    lme_te(L,1) = nanmean(TEeffM(lid));
    lme_te(L,2) = nanmean(TEeff_L(lid));
    lme_te(L,3) = nanmean(TEeff_HTL(lid));
    lme_te(L,4) = nanmean(TEeff_HTLd(lid));
    lme_te(L,5) = nanmean(TEeff_LTL(lid));
    lme_te(L,6) = nanmean(TEeff_LTLd(lid));
    
end

lme_m = NaN*ones(ni,nj);
lme_l = lme_m;
lme_htl = lme_m;
lme_ltl = lme_m;
lme_htlD = lme_m;
lme_ltlD = lme_m;
for L=1:66
    lid = find(tlme==L);

    lme_m(lid)      = lme_te(L,1);
    lme_l(lid)      = lme_te(L,2);
    lme_htl(lid)    = lme_te(L,3);
    lme_htlD(lid)   = lme_te(L,4);
    lme_ltl(lid)    = lme_te(L,5);
    lme_ltlD(lid)   = lme_te(L,6);
end

%%
save([dpath 'TEeff_Climatol_All_fish03_' cfile '.mat'],'lme_te',...
    'lme_m','lme_l','lme_htl','lme_ltl','lme_htlD','lme_ltlD','-append');

TEM = real(lme_te(:,1).^(1/2));
TEL = real(lme_te(:,2).^(1/4));
TEHTL   = real(lme_te(:,3).^(1/3));
TEHTLd  = real(lme_te(:,4).^(1/3));
TELTL   = real(lme_te(:,5));
TELTLd  = real(lme_te(:,6));

Tab=table([1:66]',lme_te(:,2),lme_te(:,3),lme_te(:,4),lme_te(:,5),...
    lme_te(:,6),TEL,TEHTL,TEHTLd,TELTL,TELTLd,...
    'VariableNames',{'LME','TEeffL','TEeffHTLb','TEeffHTLd',...
    'TEeffLTLb','TEeffLTLd','TEL','TEHTLb','TEHTLd',...
    'TELTLb','TELTLd'});
writetable(Tab,[dpath 'LME_TEeff_clim_fished_',harv,'_' cfile '.csv'],'Delimiter',',');
save([dpath 'LME_TEeff_clim_fished_',harv,'_' cfile '.mat'],'Tab');

writetable(Tab,[dpath 'LME_TEeff_clim_fished_',harv,'_' cfile '.csv'],'Delimiter',',');

%% Figures
% M
% figure(1)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,lme_m)
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0.005 0.03]);
% hcb = colorbar('h');
% ylim(hcb,[0.005 0.03])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology TEeff M')
% stamp(cfile)
% print('-dpng',[ppath 'Clim_fished_',harv,'_LME_TEeffM.png'])

% L
figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_l)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.0005 0.003]);
hcb = colorbar('h');
ylim(hcb,[0.0005 0.003])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology TEeff L')
stamp(cfile)
print('-dpng',[ppath 'Clim_fished_',harv,'_LME_TEeffL.png'])

% HTLb
figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_htl)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.001 0.01]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Climatology TEeff HTL (bent)')
stamp(cfile)
%print('-dpng',[ppath 'Clim_fished_',harv,'_LME_TEeffHTL.png'])

% HTLd
figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_htlD)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.005 0.045]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Climatology TEeff HTL (det)')
stamp(cfile)
print('-dpng',[ppath 'Clim_fished_',harv,'_LME_TEeffHTLd.png'])

% LTLb
figure(5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_ltl)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.35]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Climatology TEeff LTL (bent)')
stamp(cfile)
%print('-dpng',[ppath 'Clim_fished_',harv,'_LME_TEeffLTL.png'])

% LTLd
figure(6)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_ltlD)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.02 0.12]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Climatology TEeff LTL (det)')
stamp(cfile)
print('-dpng',[ppath 'Clim_fished_',harv,'_LME_TEeffLTLd.png'])


%% M
% figure(7)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(lme_m.^(1/2)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0 0.4]);
% hcb = colorbar('h');
% ylim(hcb,[0 0.4])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology TEeff M')
% stamp(cfile)
% print('-dpng',[ppath 'Clim_fished_',harv,'_LME_TEeffM_converted.png'])

% L
figure(8)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(lme_l.^(1/4)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.35]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Climatology TE L')
stamp(cfile)
print('-dpng',[ppath 'Clim_fished_',harv,'_LME_TEeffL_converted.png'])

% HTL
figure(9)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(lme_htl.^(1/3)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.35]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Climatology TE HTL (bent)')
stamp(cfile)
%print('-dpng',[ppath 'Clim_fished_',harv,'_LME_TEeffHTL_converted.png'])

% HTLd
figure(10)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(lme_htlD.^(1/3)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.35]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Climatology TE HTL (det)')
stamp(cfile)
print('-dpng',[ppath 'Clim_fished_',harv,'_LME_TEeffHTLd_converted.png'])

% LLT
figure(11)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(lme_ltl))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.4]);
hcb = colorbar('h');
ylim(hcb,[0.05 0.35])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology TE LTL (bent)')
stamp(cfile)
%print('-dpng',[ppath 'Clim_fished_',harv,'_LME_TEeffLTL_converted.png'])

%% LLTd
figure(12)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(lme_ltlD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.15]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Climatology TE LTL (det)')
stamp(cfile)
print('-dpng',[ppath 'Clim_fished_',harv,'_LME_TEeffLTLd_converted.png'])

