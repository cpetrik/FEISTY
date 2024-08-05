%% Calc mean nat+pred mort of F for Sophie

clear
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
exp = 'Climatol_All_fish03_';

fpath = ['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/',cfile,'/Climatology/'];

%%
% Death rates (g m-2 d-1): Sf.die 
% Predation rates (m-2 d-1): Sf.pred 
% Nat mort biom (g m-2 d-1): mort = nmort .* bio
% Natural mortality rates (m-2 d-1): Sf.nmort
% Fishing yield (g m-2 d-1): caught
% Fishing mort (m-2 d-1): fmort = caught ./ bio;

load([fpath,'Means_die_nmort_yield_',exp,cfile,'.mat'],'time','lyr',...
    'mf_die','mf_mdie','mf_tmdie','mf_mort','mf_mmort','mf_tmmort',...
    'mf_yield','mf_myield','mf_tmyield',...
    'sf_die','sf_mdie','sf_tmdie','sf_mort','sf_mmort','sf_tmmort');

%%
load([fpath,'Means_',exp,cfile,'.mat'],...
    'sf_bio','sf_mean','sf_tmean',...
    'mf_bio','mf_mean','mf_tmean');

%% Approximate biomass-specific rates
%predation
sf_pred = sf_die ./ sf_bio;
sf_smpred = sf_mdie ./ sf_mean;
sf_tmpred = sf_tmdie ./ sf_tmean;

mf_pred = mf_die ./ mf_bio;
mf_smpred = mf_mdie ./ mf_mean;
mf_tmpred = mf_tmdie ./ mf_tmean;

f_pred = sf_pred + mf_pred;
f_smpred = sf_smpred + mf_smpred;
f_tmpred = sf_tmpred + mf_tmpred;

%nat mort
sf_mort = sf_mort ./ sf_bio;
sf_smmort = sf_mmort ./ sf_mean;
sf_tmmort = sf_tmmort ./ sf_tmean;

mf_mort = mf_mort ./ mf_bio;
mf_smmort = mf_mmort ./ mf_mean;
mf_tmmort = mf_tmmort ./ mf_tmean;

f_mort = sf_mort + mf_mort;
f_smmort = sf_smmort + mf_smmort;
f_tmmort = sf_tmmort + mf_tmmort;

%fishing mort
f_yield = mf_yield ./ mf_bio;
f_smyield = mf_myield ./ mf_mean;
f_tmyield = mf_tmyield ./ mf_tmean;

% All Natural
sf_nat = sf_pred + sf_mort;
sf_smnat = sf_smpred + sf_smmort;
sf_tmnat = sf_tmpred + sf_tmmort;

mf_nat = mf_pred + mf_mort;
mf_smnat = mf_smpred + mf_smmort;
mf_tmnat = mf_tmpred + mf_tmmort;

f_nat = sf_nat + mf_nat;
f_smnat = sf_smnat + mf_smnat;
f_tmnat = sf_tmnat + mf_tmnat;

% All Natural + Fishing
mf_all = mf_pred + mf_mort + f_yield;
mf_small = mf_smpred + mf_smmort + f_smyield;
mf_tmall = mf_tmpred + mf_tmmort + f_tmyield;

%% Grid
Pdir = '/Volumes/petrik-lab/Feisty/GCM_Data/ESM26_hist/';
pp = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

ppath = [pp cfile '/Climatol/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

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

load coastlines;                     

% colors
load('MyColormaps.mat')
cm9=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...  %b
    0 0 0];...      %black
    
%% Plots in space
Bsf=NaN*ones(ni,nj);
Bmf=NaN*ones(ni,nj);
Msf=NaN*ones(ni,nj);
Mmf=NaN*ones(ni,nj);
MortF=NaN*ones(ni,nj);
AllMort=NaN*ones(ni,nj);

Bsf(ID)=sf_mean;
Bmf(ID)=mf_mean;
Msf(ID)=sf_smnat;
Mmf(ID)=mf_smnat;
MortF(ID)=f_smnat;
AllMort(ID)=mf_small;

%% Plots

%Time series
figure(1)
subplot(2,2,1)
plot(time(lyr),sf_tmnat(lyr),'k')
set(gca,'XTick',lyr(2:2:12),'XTickLabel',2:2:12)
xlabel('Month')
ylabel('SF Natural mortality (d^-^1)')

subplot(2,2,2)
plot(time(lyr),mf_tmnat(lyr),'k')
set(gca,'XTick',lyr(2:2:12),'XTickLabel',2:2:12)
xlabel('Month')
ylabel('MF Natural mortality (d^-^1)')

subplot(2,2,3)
plot(time(lyr),f_tmnat(lyr),'k')
set(gca,'XTick',lyr(2:2:12),'XTickLabel',2:2:12)
xlabel('Month')
ylabel('TotF Natural mortality (d^-^1)')

subplot(2,2,4)
plot(time(lyr),mf_tmall(lyr),'k')
set(gca,'XTick',lyr(2:2:12),'XTickLabel',2:2:12)
xlabel('Month')
ylabel('MF Nat + Fishing mortality (d^-^1)')
print('-dpng',[ppath exp 'timeseries_forage_mort.png'])


%% Map
% Mort rates
figure(2)
% SF Nat
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Msf)
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.1]);
%colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
colorbar('Position',[0.01 0.5 0.45 0.05],'orientation','horizontal')
title('SF Natural mortality (d^-^1)')

% TotF Nat
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,MortF)
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.1]);
title('TotF Natural mortality (d^-^1)')

% MF Nat
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Mmf)
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.02]);
colorbar('Position',[0.51 0.5 0.45 0.05],'orientation','horizontal')
title('MF Natural mortality (d^-^1)')

% MF All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,AllMort)
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.02]);
title('MF Nat + Fishing mortality (d^-^1)')
print('-dpng',[ppath exp 'global_forage_mort.png'])


%% Biom in mgC/m3
%gWW/m2 -> mgC/m3: (1/9)gWW->gC; 1e3:gC->mgC; (1/100)m2->m3 
Bsf(ID) = sf_mean .* (1/9) .* 1e3 .* (1/100);
Bmf(ID) = mf_mean .* (1/9) .* 1e3 .* (1/100);
AllF = Bsf+Bmf;

figure(3)
% SF 
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Bsf)
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1]);
colorbar('Position',[0.02 0.5 0.45 0.05],'orientation','horizontal')
title('mean SF (mgC m^-^3)')

% 
%subplot('Position',[0 0 0.5 0.5])

% MF 
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Bmf)
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 10]);
colorbar('Position',[0.52 0.5 0.45 0.05],'orientation','horizontal')
title('mean MF (mgC m^-^3)')

% Total
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,AllF)
cmocean('tempo')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 10]);
title('mean TotalF (mgC m^-^3)')

print('-dpng',[ppath exp 'global_forage_biom_mgCm3.png'])




