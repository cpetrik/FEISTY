% Visualize difference between
% ESM2M Hindcast of 1951-2000 and Forecast of 2051-2100
% For ensemble members that resulted in largest and smallest changes

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

harv = 'All_fish03';

%cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
lpath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([lpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_VarMaxMinDiffSims.mat']);

nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'Hist_Fore_All_fish03_ensem6_mid_kt3_maxmin_pdiffs.mat'],...
    'pnames','snames');

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

cmBP=cbrewer('seq','BuPu',50,'PCHIP');

%% max & min change side-by-side
xaid = pstats(1,1);
xfid = pstats(3,1);
xpid = pstats(5,1);
xdid = pstats(7,1);
xbid = pstats(9,1);
xapath = pnames{xaid}; xafile = snames{xaid};
xfpath = pnames{xfid}; xffile = snames{xfid};
xppath = pnames{xpid}; xpfile = snames{xpid};
xdpath = pnames{xdid}; xdfile = snames{xdid};
xbpath = pnames{xbid}; xbfile = snames{xbid};

naid = pstats(2,1);
nfid = pstats(4,1);
npid = pstats(6,1);
ndid = pstats(8,1);
nbid = pstats(10,1);
napath = pnames{naid}; nafile = snames{naid};
nfpath = pnames{nfid}; nffile = snames{nfid};
nppath = pnames{npid}; npfile = snames{npid};
ndpath = pnames{ndid}; ndfile = snames{ndid};
nbpath = pnames{nbid}; nbfile = snames{nbid};

%% MAX ----------------------------------------------------
%%% Singles
% Hindcast
load([xfpath '/Historic_All_fish03_Means_' xffile '.mat'],...
    'sf_mean50','mf_mean50');
load([xppath '/Historic_All_fish03_Means_' xpfile '.mat'],...
    'sp_mean50','mp_mean50','lp_mean50');
load([xdpath '/Historic_All_fish03_Means_' xdfile '.mat'],...
    'sd_mean50','md_mean50','ld_mean50');
load([xbpath '/Historic_All_fish03_Means_' xbfile '.mat'],...
    'b_mean50');

[ni,nj]=size(geolon_t);
xhF=NaN*ones(ni,nj);
xhP=NaN*ones(ni,nj);
xhD=NaN*ones(ni,nj);
xhB =NaN*ones(ni,nj);
xhF(grid(:,1))=sf_mean50+mf_mean50;
xhP(grid(:,1))=sp_mean50+mp_mean50+lp_mean50;
xhD(grid(:,1))=sd_mean50+md_mean50+ld_mean50;
xhB(grid(:,1))=b_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50 b_mean50

% Forecast
load([xfpath '/Forecast_All_fish03_Means_' xffile '.mat'],...
    'sf_mean50','mf_mean50');
load([xppath '/Forecast_All_fish03_Means_' xpfile '.mat'],...
    'sp_mean50','mp_mean50','lp_mean50');
load([xdpath '/Forecast_All_fish03_Means_' xdfile '.mat'],...
    'sd_mean50','md_mean50','ld_mean50');
load([xbpath '/Forecast_All_fish03_Means_' xbfile '.mat'],...
    'b_mean50');

[ni,nj]=size(geolon_t);
xcF=NaN*ones(ni,nj);
xcP=NaN*ones(ni,nj);
xcD=NaN*ones(ni,nj);
xcB =NaN*ones(ni,nj);
xcF(grid(:,1))=sf_mean50+mf_mean50;
xcP(grid(:,1))=sp_mean50+mp_mean50+lp_mean50;
xcD(grid(:,1))=sd_mean50+md_mean50+ld_mean50;
xcB(grid(:,1))=b_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50 b_mean50

xpdiffF = (xcF-xhF) ./ xhF;
xpdiffP = (xcP-xhP) ./ xhP;
xpdiffD = (xcD-xhD) ./ xhD;
xpdiffB = (xcB-xhB) ./ xhB;


%%% All
%Hindcast
load([xapath '/Historic_All_fish03_Means_' xafile '.mat'],...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50');

xhAll=NaN*ones(ni,nj);
xhAll(ID)=sf_mean50+sp_mean50+sd_mean50+mf_mean50+mp_mean50+md_mean50+...
    lp_mean50+ld_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50

load([xapath '/Forecast_All_fish03_Means_' xafile '.mat'],...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50');

xcAll=NaN*ones(ni,nj);
xcAll(ID)=sf_mean50+sp_mean50+sd_mean50+mf_mean50+mp_mean50+md_mean50+...
    lp_mean50+ld_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50

xpdiffAll = (xcAll-xhAll) ./ xhAll;

%% MIN ----------------------------------------------------
%FEISTY Output
% Hindcast
load([nfpath '/Historic_All_fish03_Means_' nffile '.mat'],...
    'sf_mean50','mf_mean50');
load([nppath '/Historic_All_fish03_Means_' npfile '.mat'],...
    'sp_mean50','mp_mean50','lp_mean50');
load([ndpath '/Historic_All_fish03_Means_' ndfile '.mat'],...
    'sd_mean50','md_mean50','ld_mean50');
load([nbpath '/Historic_All_fish03_Means_' nbfile '.mat'],...
    'b_mean50');

[ni,nj]=size(geolon_t);
nhF=NaN*ones(ni,nj);
nhP=NaN*ones(ni,nj);
nhD=NaN*ones(ni,nj);
nhB =NaN*ones(ni,nj);
nhF(grid(:,1))=sf_mean50+mf_mean50;
nhP(grid(:,1))=sp_mean50+mp_mean50+lp_mean50;
nhD(grid(:,1))=sd_mean50+md_mean50+ld_mean50;
nhB(grid(:,1))=b_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50 b_mean50

% Forecast
load([nfpath '/Forecast_All_fish03_Means_' nffile '.mat'],...
    'sf_mean50','mf_mean50');
load([nppath '/Forecast_All_fish03_Means_' npfile '.mat'],...
    'sp_mean50','mp_mean50','lp_mean50');
load([ndpath '/Forecast_All_fish03_Means_' ndfile '.mat'],...
    'sd_mean50','md_mean50','ld_mean50');
load([nbpath '/Forecast_All_fish03_Means_' nbfile '.mat'],...
    'b_mean50');

[ni,nj]=size(geolon_t);
ncF=NaN*ones(ni,nj);
ncP=NaN*ones(ni,nj);
ncD=NaN*ones(ni,nj);
ncB =NaN*ones(ni,nj);
ncF(grid(:,1))=sf_mean50+mf_mean50;
ncP(grid(:,1))=sp_mean50+mp_mean50+lp_mean50;
ncD(grid(:,1))=sd_mean50+md_mean50+ld_mean50;
ncB(grid(:,1))=b_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50 b_mean50

npdiffF = (ncF-nhF) ./ nhF;
npdiffP = (ncP-nhP) ./ nhP;
npdiffD = (ncD-nhD) ./ nhD;
npdiffB = (ncB-nhB) ./ nhB;


% All
%Hindcast
load([napath '/Historic_All_fish03_Means_' nafile '.mat'],...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50');

nhAll=NaN*ones(ni,nj);
nhAll(ID)=sf_mean50+sp_mean50+sd_mean50+mf_mean50+mp_mean50+md_mean50+...
    lp_mean50+ld_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50

load([napath '/Forecast_All_fish03_Means_' nafile '.mat'],...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50');

ncAll=NaN*ones(ni,nj);
ncAll(ID)=sf_mean50+sp_mean50+sd_mean50+mf_mean50+mp_mean50+md_mean50+...
    lp_mean50+ld_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50

npdiffAll = (ncAll-nhAll) ./ nhAll;

%% Forage
%Max
figure(5)
% Hist
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xhF))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast F');
% Fore
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xcF))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Forecast F');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('maximum Forecast - Hindcast F');
print('-dpng',[pp 'Hist_Fore_' harv '_global_maxF_3plot.png'])

%Min
figure(6)
% Hist
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(nhF))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast F');
% Fore
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(ncF))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Forecast F');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('minimum Forecast - Hindcast F');
print('-dpng',[pp 'Hist_Fore_' harv '_global_minF_3plot.png'])

%% Large Pel
%Max
figure(7)
% Hist
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xhP))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast P');
% Fore
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xcP))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Forecast P');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('maximum Forecast - Hindcast P');
print('-dpng',[pp 'Hist_Fore_' harv '_global_maxP_3plot.png'])

%Min
figure(8)
% Hist
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(nhP))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast P');
% Fore
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(ncP))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Forecast P');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('minimum Forecast - Hindcast P');
print('-dpng',[pp 'Hist_Fore_' harv '_global_minP_3plot.png'])

%% Dem
%Max
figure(9)
% Hist
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xhD))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast D');
% Fore
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xcD))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Forecast D');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('maximum Forecast - Hindcast D');
print('-dpng',[pp 'Hist_Fore_' harv '_global_maxD_3plot.png'])

%Min
figure(10)
% Hist
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(nhD))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast D');
% Fore
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(ncD))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Forecast D');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('minimum Forecast - Hindcast D');
print('-dpng',[pp 'Hist_Fore_' harv '_global_minD_3plot.png'])

%% All
%Max
figure(13)
% Hist
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xhAll))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast All');
% Fore
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xcAll))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Forecast All');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('maximum Forecast - Hindcast All');
print('-dpng',[pp 'Hist_Fore_' harv '_global_maxAll_3plot.png'])

%Min
figure(14)
% Hist
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(nhAll))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast All');
% Fore
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(ncAll))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Forecast All');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('minimum Forecast - Hindcast All');
print('-dpng',[pp 'Hist_Fore_' harv '_global_minAll_3plot.png'])

%% Bent
%Max
figure(15)
% Hist
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xhB))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast B');
% Fore
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xcB))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Forecast B');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('maximum Forecast - Hindcast B');
print('-dpng',[pp 'Hist_Fore_' harv '_global_maxB_3plot.png'])

%Min
figure(16)
% Hist
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(nhB))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast B');
% Fore
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(ncB))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Forecast B');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('minimum Forecast - Hindcast B');
print('-dpng',[pp 'Hist_Fore_' harv '_global_minB_3plot.png'])

%% Subplots of types together ----------------------------------------------------
figure(17)
subplot(5,2,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('maximum Forecast - Hindcast F');

subplot(5,2,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('minimum Forecast - Hindcast F');


subplot(5,2,3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('maximum Forecast - Hindcast P');

subplot(5,2,4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('minimum Forecast - Hindcast P');


subplot(5,2,5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('maximum Forecast - Hindcast D');

subplot(5,2,6)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('minimum Forecast - Hindcast D');


subplot(5,2,7)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('maximum Forecast - Hindcast All');

subplot(5,2,8)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('minimum Forecast - Hindcast All');


subplot(5,2,9)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('maximum Forecast - Hindcast B');

subplot(5,2,10)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('minimum Forecast - Hindcast B');

print('-dpng',[pp 'Hist_Fore_' harv '_global_maxmin_10plot.png'])

%% Just 3 functional types
figure(18)
%1
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
%title('maximum Forecast - Hindcast F');
text(-0.75,1.75,'max \DeltaF')
%text(-2.75,1.75,'A')

%2
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
%title('minimum Forecast - Hindcast F');
text(-0.75,1.75,'min \DeltaF')
%text(-2.75,1.75,'B')

%3
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
colorbar('Position',[0.8 0.25 0.025 0.5],'orientation','vertical')
%title('maximum Forecast - Hindcast P');
text(-0.75,1.75,'max \DeltaP')
%text(-2.75,1.75,'C')

%4
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
%title('minimum Forecast - Hindcast P');
text(-0.75,1.75,'min \DeltaP')
%text(-2.75,1.75,'D')

%5
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
%colorbar('Position',[0.4375 0.32 0.34 0.025],'orientation','horizontal','AxisLocation','in')
%title('maximum Forecast - Hindcast D');
text(-0.75,1.75,'max \DeltaD')
%text(-2.75,1.75,'E')

%6
subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
%title('minimum Forecast - Hindcast D');
text(-0.75,1.75,'min \DeltaD')
%text(-2.75,1.75,'F')

print('-dpng',[pp 'Hist_Fore_' harv '_global_maxmin_6plot.png'])

%% Just the diff between max & min
figure(1)
%1
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*(npdiffF-xpdiffF))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 50]);
set(gcf,'renderer','painters')
%title('maximum Forecast - Hindcast F');
text(-0.75,1.75,'min-max F')
%text(-2.75,1.75,'A')

%2
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*(npdiffP-xpdiffP))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 50]);
set(gcf,'renderer','painters')
%title('minimum Forecast - Hindcast F');
text(-0.75,1.75,'min-max P')
%text(-2.75,1.75,'B')

%3
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*(npdiffD-xpdiffD))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 50]);
set(gcf,'renderer','painters')
%title('maximum Forecast - Hindcast P');
text(-0.75,1.75,'min-max D')
%text(-2.75,1.75,'C')

%4
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*(npdiffAll-xpdiffAll))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 50]);
set(gcf,'renderer','painters')
colorbar('Position',[0.44 0.325 0.34 0.025],'orientation','horizontal')
%title('minimum Forecast - Hindcast P');
text(-0.75,1.75,'min-max All')
%text(-2.75,1.75,'D')

%5
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*(npdiffB-xpdiffB))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 50]);
set(gcf,'renderer','painters')
%colorbar('Position',[0.4375 0.32 0.34 0.025],'orientation','horizontal','AxisLocation','in')
%title('maximum Forecast - Hindcast D');
text(-0.75,1.75,'min-max B')
%text(-2.75,1.75,'E')

%6
% subplot('Position',[0.41 0.06 0.4 0.3])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,100*npdiffD)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-50 50]);
% set(gcf,'renderer','painters')
% %title('minimum Forecast - Hindcast D');
% text(-0.75,1.75,'min \DeltaD')
% %text(-2.75,1.75,'F')

print('-dpng',[pp 'Hist_Fore_' harv '_global_maxmin_5plot.png'])
