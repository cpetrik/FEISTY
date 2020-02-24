% Visualize difference between
% ESM2M Hindcast of 1951-2000 
% and Forecast of 2051-2100

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% NPP and zoop
load([bpath 'cobalt_det_temp_zoop_npp_means.mat']);

% molN/m2/s --> g/m2/d
mzloss_hist = mzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_hist = lzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mzprod_hist = mzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_hist = lzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_hist = npp_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_hist = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
ptemp_hist = ptemp_mean_hist - 273;

mzloss_fore = mzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_fore = lzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mzprod_fore = mzprod_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_fore = lzprod_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_fore = npp_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_fore = det_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
ptemp_fore = ptemp_mean_fore - 273;

zloss_hist = mzloss_hist + lzloss_hist;
zloss_fore = mzloss_fore + lzloss_fore;

l10ZlDet_hist = log10(zloss_hist./det_hist);
l10ZlDet_fore = log10(zloss_fore./det_fore);

ZlDet_hist = (zloss_hist./det_hist);
ZlDet_fore = (zloss_fore./det_fore);

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
%fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

% Hindcast
load([fpath 'Means_Historic_' harv '_' cfile '.mat'],...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50','b_mean50');

[hi,hj]=size(geolon_t);
Hsf=NaN*ones(hi,hj);
Hsp=NaN*ones(hi,hj);
Hsd=NaN*ones(hi,hj);
Hmf=NaN*ones(hi,hj);
Hmp=NaN*ones(hi,hj);
Hmd=NaN*ones(hi,hj);
Hlp=NaN*ones(hi,hj);
Hld=NaN*ones(hi,hj);
Hb =NaN*ones(hi,hj);
Hsf(grid(:,1))=sf_mean50;
Hsp(grid(:,1))=sp_mean50;
Hsd(grid(:,1))=sd_mean50;
Hmf(grid(:,1))=mf_mean50;
Hmp(grid(:,1))=mp_mean50;
Hmd(grid(:,1))=md_mean50;
Hlp(grid(:,1))=lp_mean50;
Hld(grid(:,1))=ld_mean50;
Hb(grid(:,1)) =b_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50 b_mean50

% Forecast
load([fpath 'Means_fore_',harv,'_' cfile '.mat'],...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50','b_mean50');

[ni,nj]=size(geolon_t);
Csf=NaN*ones(ni,nj);
Csp=NaN*ones(ni,nj);
Csd=NaN*ones(ni,nj);
Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);
Cb =NaN*ones(ni,nj);
Csf(ID)=sf_mean50;
Csp(ID)=sp_mean50;
Csd(ID)=sd_mean50;
Cmf(ID)=mf_mean50;
Cmp(ID)=mp_mean50;
Cmd(ID)=md_mean50;
Clp(ID)=lp_mean50;
Cld(ID)=ld_mean50;
Cb(ID) =b_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50 b_mean50

cF = Csf+Cmf;
cP = Csp+Cmp+Clp;
cD = Csd+Cmd+Cld;
cS = Csp+Csf+Csd;
cM = Cmp+Cmf+Cmd;
cL = Clp+Cld;

hF = Hsf+Hmf;
hP = Hsp+Hmp+Hlp;
hD = Hsd+Hmd+Hld;
hS = Hsp+Hsf+Hsd;
hM = Hmp+Hmf+Hmd;
hL = Hlp+Hld;

%% save hist and fore together
save([fpath 'Means_hist_fore_',harv,'_cobalt_' cfile '.mat'],'mzloss_hist',...
    'mzprod_hist','lzloss_hist','lzprod_hist','npp_hist','det_hist',...
    'ptemp_hist','mzloss_fore','mzprod_fore','lzloss_fore','lzprod_fore',...
    'npp_fore','det_fore','ptemp_fore','cF','cP','cD','cS','cM','cL','hF',...
    'hP','hD','hS','hM','hL');

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%%
cAll = cF+cP+cD;
cFracF = cF ./ cAll;
cFracP = cP ./ cAll;
cFracD = cD ./ cAll;

hAll = hF+hP+hD;
hFracF = hF ./ cAll;
hFracP = hP ./ cAll;
hFracD = hD ./ cAll;

pdiffF = (cF-hF) ./ hF;
pdiffP = (cP-hP) ./ hP;
pdiffD = (cD-hD) ./ hD;
pdiffB = (Cb-Hb) ./ Hb;
pdiffAll = (cAll-hAll) ./ hAll;

pdiffPf = (cFracP-hFracP) ./ hFracP;
pdiffFf = (cFracF-hFracF) ./ hFracF;
pdiffDf = (cFracD-hFracD) ./ hFracD;

diffPf = (cFracP-hFracP);
diffFf = (cFracF-hFracF);
diffDf = (cFracD-hFracD);

%% Maps
% pdiffs subplot 3 vert
figure(1)
%F
subplot('Position',[0.25 0.66 0.48 0.33])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,100*pdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
%colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
text(-0.1,1.5,'Forage fish','FontWeight','bold','HorizontalAlignment','center');

subplot('Position',[0.25 0.33 0.475 0.33])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,100*pdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.75 0.34 0.02 0.33],'orientation','vertical')
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
text(-0.1,1.5,'Large pelagics','FontWeight','bold','HorizontalAlignment','center');

subplot('Position',[0.25 0.0 0.475 0.33])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,100*pdiffD)
cmocean('balance')
%colormap(cmR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
%colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
text(-0.1,1.5,'Demersals','FontWeight','bold','HorizontalAlignment','center');
print('-dpng',[pp 'Hist_Fore_' harv '_global_pdiff_types_vert.png'])

%% diffs of fracs type out of total subplot 3 vert
figure(2)
subplot('Position',[0.25 0.66 0.48 0.33])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffFf)
cmocean('balance')
%colormap(cmR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.25 1.25]);
%colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
text(-0.1,1.5,'Forage fish','FontWeight','bold','HorizontalAlignment','center');

subplot('Position',[0.25 0.33 0.475 0.33])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffPf)
cmocean('balance')
%colormap(cmR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.25 1.25]);
colorbar('Position',[0.75 0.34 0.02 0.33],'orientation','vertical','Ticks',[-1 0 1])
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
text(-0.1,1.5,'Large pelagics','FontWeight','bold','HorizontalAlignment','center');

subplot('Position',[0.25 0.0 0.475 0.33])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffDf)
cmocean('balance')
%colormap(cmR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.25 1.25]);
%colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
text(-0.1,1.5,'Demersals','FontWeight','bold','HorizontalAlignment','center');
print('-dpng',[pp 'Hist_Fore_',harv,'_global_f.png'])

%% Diff Fracs Arctic and Antarctic individually
figure(3)
subplot('Position',[0 0.65 0.5 0.3])
axesm ('ortho','MapLatLimit',[60 90],'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffFf)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('\Delta  F / All');

% subplot('Position',[0 0.325 0.5 0.3])
% axesm ('ortho','MapLatLimit',[60 90],'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,diffPf)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% colorbar
% set(gcf,'renderer','painters')
% title('\Delta  P / All');

% Diff
subplot('Position',[0 0.02 0.5 0.3])
axesm ('ortho','MapLatLimit',[60 90],'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffDf)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('\Delta D / All');


subplot('Position',[0.5 0.65 0.5 0.3])
axesm ('ortho','MapLatLimit',[-90 -60],'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffFf)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('\Delta  F / All');

% subplot('Position',[0.5 0.325 0.5 0.3])
% axesm ('ortho','MapLatLimit',[-90 -60],'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,diffPf)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% colorbar
% set(gcf,'renderer','painters')
% title('\Delta  P / All');

subplot('Position',[0.5 0.02 0.5 0.3])
axesm ('ortho','MapLatLimit',[-90 -60],'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,diffDf)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('\Delta D / All');
print('-dpng',[pp 'Hist_Fore_',harv,'_arctic_fracs_noP.png'])

%% pDiff Arctic and Antarctic individually no P
figure(4)
subplot('Position',[0 0.65 0.5 0.3])
axesm ('ortho','MapLatLimit',[60 90],'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('F');

% subplot('Position',[0 0.325 0.5 0.3])
% axesm ('ortho','MapLatLimit',[60 90],'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,100*pdiffP)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% colorbar
% set(gcf,'renderer','painters')
% title('P');

% Diff
subplot('Position',[0 0.02 0.5 0.3])
axesm ('ortho','MapLatLimit',[60 90],'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,100*pdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('D');


subplot('Position',[0.5 0.65 0.5 0.3])
axesm ('ortho','MapLatLimit',[-90 -60],'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('F');

% subplot('Position',[0.5 0.325 0.5 0.3])
% axesm ('ortho','MapLatLimit',[-90 -60],'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,100*pdiffP)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% colorbar
% set(gcf,'renderer','painters')
% title('P');

subplot('Position',[0.5 0.02 0.5 0.3])
axesm ('ortho','MapLatLimit',[-90 -60],'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,100*pdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('D');
print('-dpng',[pp 'Hist_Fore_',harv,'_arctic_pdiff_noP.png'])

%% pDiff Arctic and Antarctic individually
figure(5)
subplot('Position',[0 0.67 0.5 0.275])
axesm ('ortho','MapLatLimit',[60 90],'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('F');

subplot('Position',[0 0.335 0.5 0.275])
axesm ('ortho','MapLatLimit',[60 90],'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('P');

% Diff
subplot('Position',[0 0.02 0.5 0.275])
axesm ('ortho','MapLatLimit',[60 90],'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,100*pdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('D');


subplot('Position',[0.5 0.67 0.5 0.275])
axesm ('ortho','MapLatLimit',[-90 -60],'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('F');

subplot('Position',[0.5 0.335 0.5 0.275])
axesm ('ortho','MapLatLimit',[-90 -60],'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('P');

subplot('Position',[0.5 0.02 0.5 0.275])
axesm ('ortho','MapLatLimit',[-90 -60],'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,100*pdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('D');
print('-dpng',[pp 'Hist_Fore_',harv,'_arctic_pdiff.png'])
