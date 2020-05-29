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


%% NPP and zoop
load([bpath 'cobalt_det_temp_zoop_npp_means.mat'],'det_mean_hist','det_mean_fore');

% molN/m2/s --> g/m2/d
det_hist = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_fore = det_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

%% Ensemble psets
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

nfile = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([nfile 'simnames_ensem6_mid_kt3_bestAIC_multFup_multPneg.mat'])

pdF = nan(length(ID),length(snames)+1);
pdP = pdF;
pdD = pdF;
pdA = pdF;
pdM = pdF;
pdL = pdF;
pdBb = pdF;
pdBp = pdF;

%%
for j = 1:length(snames)
    simname = snames{j};
    fname = pnames{j};
    
    %% Historic Last 50 year means
    load([fname '/Historic_All_fish03_Means_prod_' simname '.mat'],...
        'sf_prod50','sp_prod50','sd_prod50',...
        'mf_prod50','mp_prod50','md_prod50',...
        'lp_prod50','ld_prod50');
    load([fname '/Historic_All_fish03_Means_' simname '.mat'],...
        'b_mean50');
    
    Hsf=sf_prod50;
    Hsp=sp_prod50;
    Hsd=sd_prod50;
    Hsa=mf_prod50;
    Hsm=mp_prod50;
    Hsl=md_prod50;
    Hlp=lp_prod50;
    Hld=ld_prod50;
    Hbp = det_hist(ID) .* 0.075;
    Hbb = b_mean50;
    
    hF = Hsf+Hsa;
    hP = Hsp+Hsm+Hlp;
    hD = Hsd+Hsl+Hld;
    hS = Hsp+Hsf+Hsd;
    hM = Hsm+Hsa+Hsl;
    hL = Hlp+Hld;
    hBb = Hbb;
    hBp = Hbp;
    hAll = hF+hP+hD;
    
    clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50
    clear Hsf Hsp Hsd Hmf Hmp Hmd Hlp Hld Hbp Hbb
    
    %% Forecast Last 50 year means
    load([fname '/Forecast_All_fish03_Means_prod_' simname '.mat'],...
        'sf_prod50','sp_prod50','sd_prod50',...
        'mf_prod50','mp_prod50','md_prod50',...
        'lp_prod50','ld_prod50');
    load([fname '/Forecast_All_fish03_Means_' simname '.mat'],...
        'b_mean50');
    
    Csf=sf_prod50;
    Csp=sp_prod50;
    Csd=sd_prod50;
    Cmf=mf_prod50;
    Cmp=mp_prod50;
    Cmd=md_prod50;
    Clp=lp_prod50;
    Cld=ld_prod50;
    % calculate benthic production = detritus flux * benthic efficiency
    Cbp = det_hist(ID) .* 0.075;
    Cbb = b_mean50;
    
    cF = Csf+Cmf;
    cP = Csp+Cmp+Clp;
    cD = Csd+Cmd+Cld;
    cS = Csp+Csf+Csd;
    cM = Cmp+Cmf+Cmd;
    cL = Clp+Cld;
    cBb = Cbb;
    cBp = Cbp;
    cAll = cF+cP+cD;
    
    clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50
    clear Csf Csp Csd Cmf Cmp Cmd Clp Cld Cbb Cbp
    
    %% Group percent diffs
    pdF(:,j) = (cF-hF) ./ hF;
    pdP(:,j) = (cP-hP) ./ hP;
    pdD(:,j) = (cD-hD) ./ hD;
    pdA(:,j) = (cAll-hAll) ./ hAll;
    pdM(:,j) = (cM-hM) ./ hM;
    pdL(:,j) = (cL-hL) ./ hL;
    pdBb(:,j) = (cBb-hBb) ./ hBb;
    pdBp(:,j) = (cBp-hBp) ./ hBp;
    
end

%% FEISTY Output orig pset
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
dpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

%%
load([dpath 'Prod_pdiffs_hist_fore_',harv,'_' cfile '.mat'],'pdiffL','pdiffM',...
    'pdiffF','pdiffP','pdiffD','pdiffB','pdiffAll')
% load([fpath 'Prod_pdiffs_hist_fore_',harv,'_' cfile '.mat'],'pdiffL','pdiffM',...
%     'pdiffF','pdiffP','pdiffD','pdiffB','pdiffAll')

% Hindcast
load([fpath 'Means_Historic_' harv '_prod_' cfile '.mat'],'b_mean50');
Hb = b_mean50;
clear b_mean50

% Forecast
load([fpath 'Means_fore_',harv,'_' cfile '.mat'],'b_mean50');
Cb = b_mean50;
clear b_mean50

pdiffBb = (Cb-Hb) ./ Hb;

pdF(:,j+1) = pdiffF(ID);
pdP(:,j+1) = pdiffP(ID);
pdD(:,j+1) = pdiffD(ID);
pdA(:,j+1) = pdiffAll(ID);
pdM(:,j+1) = pdiffM(ID);
pdL(:,j+1) = pdiffL(ID);
pdBb(:,j+1) = pdiffBb;
pdBp(:,j+1) = pdiffB(ID);

%% calc std and coeff var

sF = std(100*pdF,0,2);
sP = std(100*pdP,0,2);
sD = std(100*pdD,0,2);
sA = std(100*pdA,0,2);
sM = std(100*pdM,0,2);
sL = std(100*pdL,0,2);
sBb = std(100*pdBb,0,2);
sBp = std(100*pdBp,0,2);

cvF = sF ./ mean(100*pdF,2);
cvP = sP ./ mean(100*pdP,2);
cvD = sD ./ mean(100*pdD,2);
cvA = sA ./ mean(100*pdA,2);
cvM = sM ./ mean(100*pdM,2);
cvL = sL ./ mean(100*pdL,2);
cvBb = sBb ./ mean(100*pdBb,2);
cvBp = sBp ./ mean(100*pdBp,2);

%% find regions with opposite signs
% how to do without a loop?
% find neg and pos and sum

negF = (pdF<0);
negP = (pdP<0);
negD = (pdD<0);
negA = (pdA<0);
negM = (pdM<0);
negL = (pdL<0);
negB = (pdBb<0);

nF = sum(negF,2);
nP = sum(negP,2);
nD = sum(negD,2);
nA = sum(negA,2);
nM = sum(negM,2);
nL = sum(negL,2);
nB = sum(negB,2);

posF = (pdF>0);
posP = (pdP>0);
posD = (pdD>0);
posA = (pdA>0);
posM = (pdM>0);
posL = (pdL>0);
posB = (pdBb>0);

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
epath = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
efile = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
save([epath 'Prod_diff_50yr_ensem6_mid_kt3_bestAIC_multFup_multPneg.mat'],...
    'pdF','pdP','pdD','pdA','pdM','pdL','pdBp','pdBb',...
    '-append');
save([efile 'Prod_diff_50yr_ensem6_mid_kt3_bestAIC_multFup_multPneg.mat'],...
    'pdF','pdP','pdD','pdA','pdM','pdL','pdBp','pdBb',...
    '-append');

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
Agf(grid(:,1))=agF/44;
Agp(grid(:,1))=agP/44;
Agd(grid(:,1))=agD/44;
Aga(grid(:,1))=agA/44;
Agm(grid(:,1))=agM/44;
Agl(grid(:,1))=agL/44;
Agb(grid(:,1))=agB/44;

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
print('-dpng',[ppath 'Hist_Fore_',harv,'_global_pdiff_prod_std_8plot.png'])

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
print('-dpng',[ppath 'Hist_Fore_',harv,'_global_pdiff_prod_agree_8plot.png'])


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
print('-dpng',[ppath 'Hist_Fore_',harv,'_global_pdiff_prod_std_agree_8plot.png'])