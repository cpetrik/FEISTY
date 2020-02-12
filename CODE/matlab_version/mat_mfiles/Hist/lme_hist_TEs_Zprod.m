% Plot effective TEs at LME scale
% Historic
% 1860-2005, last 50 years
% Saved as mat files
% Use Zprod instead of loss

clear all
close all

dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm2m_hist/';
cdir = '/Volumes/GFDL/GCM_DATA/ESM2M_hist/';
load([gpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t','AREA_OCN');
load([gpath 'lme_mask_esm2m.mat']);
grid = csvread([gpath 'grid_csv.csv']);
load([cpath 'cobalt_det_biom_means.mat']);
load([cpath 'cobalt_npp_means.mat']);
load([cpath 'cobalt_temp_means.mat']);
load([cpath 'cobalt_zoop_biom_means.mat']); 
%lme areas?

% plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
[ni,nj]=size(geolon_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
% ENTER -100 TO MAP ORIGIN LONG

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

ID = grid(:,1);

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

%%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'TEeffDetZprod_Historic_All_fish03_' cfile '.mat']);

%% Calc LMEs
tlme = lme_mask_esm2m';

lme_te = NaN*ones(66,2);
for L=1:66
    lid = find(tlme==L);
    %TEeff
    lme_te(L,1) = nanmean(TEeffM(lid));
    lme_te(L,2) = nanmean(TEeff_L(lid));
    lme_te(L,4) = nanmean(TEeff_HTLd(lid));
    lme_te(L,6) = nanmean(TEeff_LTLd(lid));
    
end

lme_m = NaN*ones(ni,nj);
lme_l = lme_m;
lme_htlD = lme_m;
lme_ltlD = lme_m;
for L=1:66
    lid = find(tlme==L);

    lme_m(lid)      = lme_te(L,1);
    lme_l(lid)      = lme_te(L,2);
    lme_htlD(lid)   = lme_te(L,4);
    lme_ltlD(lid)   = lme_te(L,6);
end

%%
save([dpath 'TEeffDet_Historic_All_fish03_' cfile '.mat'],'lme_te',...
    'lme_m','lme_l','lme_htlD','lme_ltlD','-append');

TEM = real(lme_te(:,1).^(1/2));
TEL = real(lme_te(:,2).^(1/4));
TEHTLd  = real(lme_te(:,4).^(1/3));
TELTLd  = real(lme_te(:,6));

Tab=table([1:66]',lme_te(:,2),lme_te(:,4),lme_te(:,6),...
    TEL,TEHTLd,TELTLd,...
    'VariableNames',{'LME','TEeffL','TEeffHTLd','TEeffLTLd',...
    'TEL','TEHTLd','TELTLd'});
writetable(Tab,[dpath 'LME_TEeff_hist_',harv,'_' cfile '.csv'],'Delimiter',',');
save([dpath 'LME_TEeff_hist_',harv,'_' cfile '.mat'],'Tab');

writetable(Tab,[dpath 'LME_TEeff_hist_',harv,'_' cfile '.csv'],'Delimiter',',');

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
% title('Historic TEeff M')
% stamp(cfile)
% print('-dpng',[ppath 'Hist_fished_',harv,'_LME_TEeffM.png'])

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
title('Historic TEeff L')
stamp(cfile)
%print('-dpng',[ppath 'Hist_fished_',harv,'_LME_TEeffL.png'])

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
title('Historic TEeff HTL (det)')
stamp(cfile)
%print('-dpng',[ppath 'Hist_fished_',harv,'_LME_TEeffHTLd.png'])

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
title('Historic TEeff LTL (det)')
stamp(cfile)
%print('-dpng',[ppath 'Hist_fished_',harv,'_LME_TEeffLTLd.png'])


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
% title('Historic TEeff M')
% stamp(cfile)
% print('-dpng',[ppath 'Hist_fished_',harv,'_LME_TEeffM_converted.png'])

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
title('Historic TE L')
stamp(cfile)
%print('-dpng',[ppath 'Hist_fished_',harv,'_LME_TEeffL_converted.png'])

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
title('Historic TE HTL (det)')
stamp(cfile)
%print('-dpng',[ppath 'Hist_fished_',harv,'_LME_TEeffHTLd_converted.png'])

% LLTd
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
title('Historic TE LTL (det)')
stamp(cfile)
%print('-dpng',[ppath 'Hist_fished_',harv,'_LME_TEeffLTLd_converted.png'])

