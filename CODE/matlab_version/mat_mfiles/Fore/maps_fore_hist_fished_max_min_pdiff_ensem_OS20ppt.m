% Visualize difference between
% ESM2M Hindcast of 1951-2000 and Forecast of 2051-2100
% For ensemble members that resulted in largest and smallest changes

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'Hist_Fore_All_fish03_ensem6_mid_kt3_maxmin_pdiffs.mat'])

harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

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

% colors
cmG=cbrewer('seq','Greens',50,'PCHIP');
cmB=cbrewer('seq','Blues',50,'PCHIP');
cmP=cbrewer('seq','Purples',50,'PCHIP');
cmR=cbrewer('seq','Reds',50,'PCHIP');
cmBW=cbrewer('seq','Greys',50,'PCHIP');

%% max & min change side-by-side
xmid = big2(1,1);
xlid = big2(2,1);
xfid = big2(3,1);
xpid = big2(4,1);
xdid = big2(5,1);
xbid = big2(8,1);
xaid = big2(7,1);
xeid = big2(6,1);
xmpath = pnames{xmid}; xmfile = snames{xmid};
xlpath = pnames{xlid}; xlfile = snames{xlid};
xfpath = pnames{xfid}; xffile = snames{xfid};
xppath = pnames{xpid}; xpfile = snames{xpid};
xdpath = pnames{xdid}; xdfile = snames{xdid};
xbpath = pnames{xbid}; xbfile = snames{xbid};
xapath = pnames{xaid}; xafile = snames{xaid};
xepath = pnames{xeid}; xefile = snames{xeid};

nmid = sml2(1,1);
nlid = sml2(2,1);
nfid = sml2(3,1);
npid = sml2(4,1);
ndid = sml2(5,1);
nbid = sml2(8,1);
naid = sml2(7,1);
neid = sml2(6,1);
nmpath = pnames{nmid}; nmfile = snames{nmid};
nlpath = pnames{nlid}; nlfile = snames{nlid};
nfpath = pnames{nfid}; nffile = snames{nfid};
nppath = pnames{npid}; npfile = snames{npid};
ndpath = pnames{ndid}; ndfile = snames{ndid};
nbpath = pnames{nbid}; nbfile = snames{nbid};
napath = pnames{naid}; nafile = snames{naid};
nepath = pnames{neid}; nefile = snames{neid};

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


%%% All Pelagic
%Hindcast
load([xepath '/Historic_All_fish03_Means_' xefile '.mat'],...
    'sf_mean50','sp_mean50','mf_mean50','mp_mean50','lp_mean50');

xhPel=NaN*ones(ni,nj);
xhPel(ID)=sf_mean50+sp_mean50+mf_mean50+mp_mean50+lp_mean50;

clear sf_mean50 sp_mean50 mf_mean50 mp_mean50 lp_mean50

load([xepath '/Forecast_All_fish03_Means_' xefile '.mat'],...
    'sf_mean50','sp_mean50','mf_mean50','mp_mean50','lp_mean50');

xcPel=NaN*ones(ni,nj);
xcPel(ID)=sf_mean50+sp_mean50+mf_mean50+mp_mean50+lp_mean50;

clear sf_mean50 sp_mean50 mf_mean50 mp_mean50 lp_mean50

xpdiffPel = (xcPel-xhPel) ./ xhPel;


%%% Medium
%Hindcast
load([xmpath '/Historic_All_fish03_Means_' xmfile '.mat'],...
    'mf_mean50','mp_mean50','md_mean50');

xhM=NaN*ones(ni,nj);
xhM(ID)=mf_mean50+mp_mean50+md_mean50;

clear mf_mean50 mp_mean50 md_mean50

load([xmpath '/Forecast_All_fish03_Means_' xmfile '.mat'],...
    'mf_mean50','mp_mean50','md_mean50');

xcM=NaN*ones(ni,nj);
xcM(ID)=mf_mean50+mp_mean50+md_mean50;

clear mf_mean50 mp_mean50 md_mean50

xpdiffM = (xcM-xhM) ./ xhM;


%%% Large
%Hindcast
load([xlpath '/Historic_All_fish03_Means_' xlfile '.mat'],...
    'lp_mean50','ld_mean50');

xhL=NaN*ones(ni,nj);
xhL(ID)=lp_mean50+ld_mean50;

clear lp_mean50 ld_mean50

load([xlpath '/Forecast_All_fish03_Means_' xlfile '.mat'],...
    'lp_mean50','ld_mean50');

xcL=NaN*ones(ni,nj);
xcL(ID)=lp_mean50+ld_mean50;

clear lp_mean50 ld_mean50

xpdiffL = (xcL-xhL) ./ xhL;

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


%%% All Pelagic
%Hindcast
load([nepath '/Historic_All_fish03_Means_' nefile '.mat'],...
    'sf_mean50','sp_mean50','mf_mean50','mp_mean50','lp_mean50');

nhPel=NaN*ones(ni,nj);
nhPel(ID)=sf_mean50+sp_mean50+mf_mean50+mp_mean50+lp_mean50;

clear sf_mean50 sp_mean50 mf_mean50 mp_mean50 lp_mean50

load([nepath '/Forecast_All_fish03_Means_' nefile '.mat'],...
    'sf_mean50','sp_mean50','mf_mean50','mp_mean50','lp_mean50');

ncPel=NaN*ones(ni,nj);
ncPel(ID)=sf_mean50+sp_mean50+mf_mean50+mp_mean50+lp_mean50;

clear sf_mean50 sp_mean50 mf_mean50 mp_mean50 lp_mean50

npdiffPel = (ncPel-nhPel) ./ nhPel;


%%% Medium
%Hindcast
load([nmpath '/Historic_All_fish03_Means_' nmfile '.mat'],...
    'mf_mean50','mp_mean50','md_mean50');

nhM=NaN*ones(ni,nj);
nhM(ID)=mf_mean50+mp_mean50+md_mean50;

clear mf_mean50 mp_mean50 md_mean50

load([nmpath '/Forecast_All_fish03_Means_' nmfile '.mat'],...
    'mf_mean50','mp_mean50','md_mean50');

ncM=NaN*ones(ni,nj);
ncM(ID)=mf_mean50+mp_mean50+md_mean50;

clear mf_mean50 mp_mean50 md_mean50

npdiffM = (ncM-nhM) ./ nhM;


%%% Large
%Hindcast
load([nlpath '/Historic_All_fish03_Means_' nlfile '.mat'],...
    'lp_mean50','ld_mean50');

nhL=NaN*ones(ni,nj);
nhL(ID)=lp_mean50+ld_mean50;

clear lp_mean50 ld_mean50

load([nlpath '/Forecast_All_fish03_Means_' nlfile '.mat'],...
    'lp_mean50','ld_mean50');

ncL=NaN*ones(ni,nj);
ncL(ID)=lp_mean50+ld_mean50;

clear lp_mean50 ld_mean50

npdiffL = (ncL-nhL) ./ nhL;


%% Med
%Distr
figure(1)
% Max
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xcM))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass max \DeltaM');
% Min
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(ncM))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass min \DeltaM');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(log10(xcM)-log10(ncM)))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('difference max-min \DeltaM');
print('-dpng',[pp 'Hist_Fore_' harv '_global_distr_maxminM_3plot.png'])

%% Pdiff M
figure(2)
% Max
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*xpdiffM)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('%change max \DeltaM');
% Min
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*npdiffM)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('%change min \DeltaM');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*(xpdiffM-npdiffM))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('difference max-min \DeltaM');
print('-dpng',[pp 'Hist_Fore_' harv '_global_pdiff_maxminM_3plot.png'])

%% Large
%Distr
figure(3)
% Max
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xcL))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass max \DeltaL');
% Min
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(ncL))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass min \DeltaL');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(log10(xcL)-log10(ncL)))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('difference max-min \DeltaL');
print('-dpng',[pp 'Hist_Fore_' harv '_global_distr_maxminL_3plot.png'])

%% Forage
%Distr
figure(4)
% Max
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xcF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass max \DeltaF');
% Min
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(ncF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass min \DeltaF');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(log10(xcF)-log10(ncF)))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('difference max-min \DeltaF');
print('-dpng',[pp 'Hist_Fore_' harv '_global_distr_maxminF_3plot.png'])

%% All Pel
%Distr
figure(5)
% Max
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xcPel))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass max \DeltaPel');
% Min
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(ncPel))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass min \DeltaPel');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(log10(xcPel)-log10(ncPel)))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('difference max-min \DeltaPel');
print('-dpng',[pp 'Hist_Fore_' harv '_global_distr_maxminPel_3plot.png'])

%% Dem
%Distr
figure(6)
% Max
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xcD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass max \DeltaD');
% Min
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(ncD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass min \DeltaD');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(log10(xcD)-log10(ncD)))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('difference max-min \DeltaD');
print('-dpng',[pp 'Hist_Fore_' harv '_global_distr_maxminD_3plot.png'])

%% Large Pel
%Distr
figure(7)
% Max
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xcP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass max \DeltaP');
% Min
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(ncP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass min \DeltaP');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(log10(xcP)-log10(ncP)))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('difference max-min \DeltaP');
print('-dpng',[pp 'Hist_Fore_' harv '_global_distr_maxminP_3plot.png'])

%% All
%Distr
figure(8)
% Max
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xcAll))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass max \DeltaAll');
% Min
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(ncAll))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass min \DeltaAll');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(log10(xcAll)-log10(ncAll)))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('difference max-min \DeltaAll');
print('-dpng',[pp 'Hist_Fore_' harv '_global_distr_maxminAll_3plot.png'])

%% Bent
%Distr
figure(9)
% Max
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(xcB))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass max \DeltaB');
% Min
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(ncB))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 biomass min \DeltaB');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,(log10(xcB)-log10(ncB)))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('difference max-min \DeltaB');
print('-dpng',[pp 'Hist_Fore_' harv '_global_distr_maxminB_3plot.png'])

%% Subplots of types together ----------------------------------------------------
% Just 3 functional types
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
colorbar('Position',[0.8 0.3 0.02 0.5],'AxisLocation','out')
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

%% Just 2 size + bent
figure(19)
%1
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffM)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
%title('maximum Forecast - Hindcast F');
text(-0.75,1.75,'max \DeltaM')
%text(-2.75,1.75,'A')
%2
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffM)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
%title('minimum Forecast - Hindcast F');
text(-0.75,1.75,'min \DeltaM')
%text(-2.75,1.75,'B')
%3
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
%title('maximum Forecast - Hindcast P');
text(-0.75,1.75,'max \DeltaL')
%text(-2.75,1.75,'C')
%4
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
%title('minimum Forecast - Hindcast P');
text(-0.75,1.75,'min \DeltaL')
%text(-2.75,1.75,'D')
%5
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*xpdiffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
colorbar('Position',[0.8 0.3 0.02 0.5],'AxisLocation','out')
%title('maximum Forecast - Hindcast D');
text(-0.75,1.75,'max \DeltaB')
%text(-2.75,1.75,'E')
%6
subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*npdiffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
%title('minimum Forecast - Hindcast D');
text(-0.75,1.75,'min \DeltaB')
%text(-2.75,1.75,'F')

print('-dpng',[pp 'Hist_Fore_' harv '_global_maxmin_size_6plot.png'])

%% Pdiff All
figure(20)
% Max
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*xpdiffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
%colorbar
set(gcf,'renderer','painters')
title('%change max \DeltaAll');
% Min
subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*npdiffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
%colorbar
set(gcf,'renderer','painters')
title('%change min \DeltaAll');
% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,100*(xpdiffAll-npdiffAll))
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar
set(gcf,'renderer','painters')
title('difference max-min \DeltaAll');
print('-dpng',[pp 'Hist_Fore_' harv '_global_pdiff_maxminAll_3plot.png'])


