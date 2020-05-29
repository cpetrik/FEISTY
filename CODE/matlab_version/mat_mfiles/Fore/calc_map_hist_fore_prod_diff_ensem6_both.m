% Calc changes in production to make maps of variance across ensemble
% Ensemble 6, temp3 and orig
% recalculate benthic production = detritus flux * benthic efficiency

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

harv = 'All_fish03';

%% NPP and zoop
load([bpath 'cobalt_det_temp_zoop_npp_means.mat'],'det_mean_hist','det_mean_fore');

% molN/m2/s --> g/m2/d
det_hist = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_fore = det_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

%% Ensemble psets
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

epath = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([epath 'Prod_diff_50yr_ensem6_mid_kt3_bestAIC_multFup_multPneg.mat'],...
    'pdF','pdP','pdD','pdA','pdM','pdL','pdBp','pdBb');

dF = pdF;
dP = pdP;
dD = pdD;
dA = pdA;
dM = pdM;
dL = pdL;
dB = pdBb;

clear pdF pdP pdD pdA pdM pdL pdBp pdBb

nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'Prod_diff_50yr_ensem6_mid_samek_bestAIC_Fupneg_mult8_Pneg2_mult3.mat'],...
    'pdF','pdP','pdD','pdA','pdM','pdL','pdBp','pdBb');

dF(:,45:59) = pdF;
dP(:,45:59) = pdP;
dD(:,45:59) = pdD;
dA(:,45:59) = pdA;
dM(:,45:59) = pdM;
dL(:,45:59) = pdL;
dB(:,45:59) = pdBb;

clear pdF pdP pdD pdA pdM pdL pdBp pdBb

%% calc std and coeff var

sF = std(100*dF,0,2);
sP = std(100*dP,0,2);
sD = std(100*dD,0,2);
sA = std(100*dA,0,2);
sM = std(100*dM,0,2);
sL = std(100*dL,0,2);
sBb = std(100*dB,0,2);

cvF = sF ./ mean(100*dF,2);
cvP = sP ./ mean(100*dP,2);
cvD = sD ./ mean(100*dD,2);
cvA = sA ./ mean(100*dA,2);
cvM = sM ./ mean(100*dM,2);
cvL = sL ./ mean(100*dL,2);
cvBb = sBb ./ mean(100*dB,2);

%% find regions with opposite signs
% how to do without a loop?
% find neg and pos and sum

negF = (dF<0);
negP = (dP<0);
negD = (dD<0);
negA = (dA<0);
negM = (dM<0);
negL = (dL<0);
negB = (dB<0);

nF = sum(negF,2);
nP = sum(negP,2);
nD = sum(negD,2);
nA = sum(negA,2);
nM = sum(negM,2);
nL = sum(negL,2);
nB = sum(negB,2);

posF = (dF>0);
posP = (dP>0);
posD = (dD>0);
posA = (dA>0);
posM = (dM>0);
posL = (dL>0);
posB = (dB>0);

pF = sum(posF,2);
pP = sum(posP,2);
pD = sum(posD,2);
pA = sum(posA,2);
pM = sum(posM,2);
pL = sum(posL,2);
pB = sum(posB,2);

%agreement = max
agF = max(nF,pF);
agP = max(nP,pP);
agD = max(nD,pD);
agA = max(nA,pA);
agM = max(nM,pM);
agL = max(nL,pL);
agB = max(nB,pB);

%% save
mpath = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/';
mfile = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/';
save([mpath 'Prod_diff_50yr_ensem6_mid_both_bestAIC_multFup_multPneg.mat'],...
    'dF','dP','dD','dA','dM','dL','dB');
save([mfile 'Prod_diff_50yr_ensem6_mid_both_bestAIC_multFup_multPneg.mat'],...
    'dF','dP','dD','dA','dM','dL','dB');

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

cmBP=cbrewer('seq','BuPu',50,'PCHIP');

[hi,hj]=size(geolon_t);
Hsf=NaN*ones(hi,hj);
Hsp=NaN*ones(hi,hj);
Hsd=NaN*ones(hi,hj);
Hsa=NaN*ones(hi,hj);
Hsm=NaN*ones(hi,hj);
Hsl=NaN*ones(hi,hj);
Hsb=NaN*ones(hi,hj);
Hsf(grid(:,1))=sF;
Hsp(grid(:,1))=sP;
Hsd(grid(:,1))=sD;
Hsa(grid(:,1))=sA;
Hsm(grid(:,1))=sM;
Hsl(grid(:,1))=sL;
Hsb(grid(:,1))=sBb;

Agf=NaN*ones(hi,hj);
Agp=NaN*ones(hi,hj);
Agd=NaN*ones(hi,hj);
Aga=NaN*ones(hi,hj);
Agm=NaN*ones(hi,hj);
Agl=NaN*ones(hi,hj);
Agb=NaN*ones(hi,hj);
Agf(grid(:,1))=agF/59;
Agp(grid(:,1))=agP/59;
Agd(grid(:,1))=agD/59;
Aga(grid(:,1))=agA/59;
Agm(grid(:,1))=agM/59;
Agl(grid(:,1))=agL/59;
Agb(grid(:,1))=agB/59;

%% map std
f1 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - npp
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Hsm)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Medium','HorizontalAlignment','center')
text(-2.5,1.75,'A')

% - large
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Hsl)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Large','HorizontalAlignment','center')
text(-2.5,1.75,'B')

% - benthos
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Hsb)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Benthos','HorizontalAlignment','center')
text(-2.5,1.75,'C')

% % - benthos
% subplot('Position',[0.025 0.0 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,100*pdiffB)
% colormap(cmBP)
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'Benthos','HorizontalAlignment','center')
% text(-2.5,1.75,'D')

% - forage
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Hsf)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Forage','HorizontalAlignment','center')
text(-2.5,1.75,'E')

% - P
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Hsp)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 100]);
colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Large pelagic','HorizontalAlignment','center')
text(-2.5,1.75,'F')

% D
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Hsd)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')
text(-2.5,1.75,'G')

%- All
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Hsa)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 100]);
set(gcf,'renderer','painters')
text(0,1.75,'All','HorizontalAlignment','center')
text(-2.5,1.75,'H')
print('-dpng',[ppath 'Hist_Fore_',harv,'_global_pdiff_prod_both_std_8plot.png'])

%% map agreement
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f2.Units = 'inches';

%A - medium
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*Agm)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([50 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Medium','HorizontalAlignment','center')
text(-2.5,1.75,'A')

% - large
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*Agl)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([50 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Large','HorizontalAlignment','center')
text(-2.5,1.75,'B')

% - benthos
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*Agb)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([50 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Benthos','HorizontalAlignment','center')
text(-2.5,1.75,'C')

% % - benthos
% subplot('Position',[0.025 0.0 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,100*pdiffB)
% colormap(cmBP)
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([50 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'Benthos','HorizontalAlignment','center')
% text(-2.5,1.75,'D')

% - forage
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*Agf)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([50 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Forage','HorizontalAlignment','center')
text(-2.5,1.75,'E')

% - P
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*Agp)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([50 100]);
colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Large pelagic','HorizontalAlignment','center')
text(-2.5,1.75,'F')

% D
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*Agd)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([50 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')
text(-2.5,1.75,'G')

%- All
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*Aga)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([50 100]);
set(gcf,'renderer','painters')
text(0,1.75,'All','HorizontalAlignment','center')
text(-2.5,1.75,'H')
print('-dpng',[ppath 'Hist_Fore_',harv,'_global_pdiff_prod_both_agree_8plot.png'])


%% map std & agree together
f3 = figure('Units','inches','Position',[1 3 6.5 8]);
%f3.Units = 'inches';

%A - F
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Hsf)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 50]);
set(gcf,'renderer','painters')
text(0,1.75,{'std';'Forage'},'HorizontalAlignment','center')
text(-2.5,1.75,'A')

% - P
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Hsp)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 50]);
colorbar('Position',[0.43 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Large pelagic','HorizontalAlignment','center')
text(-2.5,1.75,'B')

% - D
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Hsd)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 50]);
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')
text(-2.5,1.75,'C')

% - All
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Hsa)
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 50]);
set(gcf,'renderer','painters')
text(0,1.75,'All','HorizontalAlignment','center')
text(-2.5,1.75,'D')

% - forage
ax1 = subplot('Position',[0.51 0.75 0.4 0.25]);
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*Agf)
colormap(ax1,flipud(cmBP))
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([50 100]);
set(gcf,'renderer','painters')
text(0,1.75,{'agree';'Forage'},'HorizontalAlignment','center')
text(-2.5,1.75,'E')

% - P
ax2 = subplot('Position',[0.51 0.5 0.4 0.25]);
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*Agp)
colormap(ax2,flipud(cmBP))
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([50 100]);
colorbar('Position',[0.925 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Large pelagic','HorizontalAlignment','center')
text(-2.5,1.75,'F')

% D
ax3 = subplot('Position',[0.51 0.25 0.4 0.25]);
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*Agd)
colormap(ax3,flipud(cmBP))
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([50 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')
text(-2.5,1.75,'G')

%- All
ax4 = subplot('Position',[0.51 0.0 0.4 0.25]);
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*Aga)
colormap(ax4, flipud(cmBP))
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([50 100]);
set(gcf,'renderer','painters')
text(0,1.75,'All','HorizontalAlignment','center')
text(-2.5,1.75,'H')
print('-dpng',[ppath 'Hist_Fore_',harv,'_global_pdiff_prod_both_std_agree_8plot.png'])





