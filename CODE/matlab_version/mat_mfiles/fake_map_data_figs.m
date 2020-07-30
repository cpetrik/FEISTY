% Fake map data for MAPP schematic

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/Historic_ESM2M/'];
%fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

ppath = '/Users/cpetrik/Dropbox/ImptDocs/Funding/NOAA/ClimatePrograms/NCAR/schematic/';

load([fpath 'Means_Historic_' harv '_' cfile '.mat'],'y',...
    'mz_tmfrac','mz_mfrac50','mz_mfrac5','mz_mfrac90','mz_mfrac','mz_ttf',...
    'mz_mtf50','mz_mtf5','mz_mtf90','mz_mtf',...
    'lz_tmfrac','lz_mfrac50','lz_mfrac5','lz_mfrac90','lz_mfrac','lz_ttf',...
    'lz_mtf50','lz_mtf5','lz_mtf90','lz_mtf');

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

[ni,nj]=size(geolon_t);
ID = grid(:,1);
nx = length(ID);

%% Pick which time period mean
% 1951-2000
%Mean fraction
Hmz_smfrac = mz_mfrac50;
Hlz_smfrac = lz_mfrac50;

%Total times it happens over time
Hmz_ttover = mz_ttf/nx;
Hlz_ttover = lz_ttf/nx;

%Total times it happens in 50 yrs * 12 mos in space
Hmz_stover = mz_mtf50 ./ (50*12);
Hlz_stover = lz_mtf50 ./ (50*12);


%% Plots in space

HFmz=NaN*ones(ni,nj);
HOmz=NaN*ones(ni,nj);
HFlz=NaN*ones(ni,nj);
HOlz=NaN*ones(ni,nj);

HFmz(ID)=(Hmz_smfrac);
HFlz(ID)=(Hlz_smfrac);

HOmz(ID)=(Hmz_stover);
HOlz(ID)=(Hlz_stover);

HFmz=rescale(HFmz,-1,1);
HOmz=rescale(HOmz,-1,1);
HFlz=rescale(HFlz,-1,1);
HOlz=rescale(HOlz,-1,1);

t1 = HFmz + HOmz;
t2 = HFmz + HOlz;
t3 = HFlz + HOmz;
t4 = HFlz + HOlz;
t5 = HOmz + HOlz;

t1=rescale(t1,-1,1);
t2=rescale(t2,-1,1);
t3=rescale(t3,-1,1);
t4=rescale(t4,-1,1);

%% fish
load([fpath 'Means_Historic_' harv '_prod_' cfile '.mat']);

% Pick which time period mean
% 1951-2000
sp_smean=sp_prod50;
sf_smean=sf_prod50;
sd_smean=sd_prod50;
mp_smean=mp_prod50;
mf_smean=mf_prod50;
md_smean=md_prod50;
lp_smean=lp_prod50;
ld_smean=ld_prod50;
b_smean=b_mean50;

Psf=NaN*ones(ni,nj);
Psp=NaN*ones(ni,nj);
Psd=NaN*ones(ni,nj);
Pmf=NaN*ones(ni,nj);
Pmp=NaN*ones(ni,nj);
Pmd=NaN*ones(ni,nj);
Plp=NaN*ones(ni,nj);
Pld=NaN*ones(ni,nj);
Pb=NaN*ones(ni,nj);

Psf(ID)=sf_prod50;
Psp(ID)=sp_prod50;
Psd(ID)=sd_prod50;
Pmf(ID)=mf_prod50;
Pmp(ID)=mp_prod50;
Pmd(ID)=md_prod50;
Plp(ID)=lp_prod50;
Pld(ID)=ld_prod50;
Pb(ID)=b_mean50;

Psf(Psf(:)<0) = nan; 
Psp(Psp(:)<0) = nan;
Psd(Psd(:)<0) = nan;
Pmf(Pmf(:)<0) = nan;
Pmp(Pmp(:)<0) = nan;
Pmd(Pmd(:)<0) = nan;
Plp(Plp(:)<0) = nan;
Pld(Pld(:)<0) = nan;
Pb(Pb(:)<0) = nan;

All = Psp+Psf+Psd+Pmp+Pmf+Pmd+Plp+Pld;
AllM = Pmp+Pmf+Pmd;
AllL = Plp+Pld;

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

cmBP=cbrewer('seq','BuPu',50,'PCHIP');
cmBP2 = cmBP;
cmBP2(11,:) = [0 0 0];

%% 4 plot of all maps
figure(1)
%1 - m frac
subplot('Position',[0.01 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,HFmz*-1)
%colormap(cmBP2)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
%colorbar('Position',[0.385 0.695 0.025 0.275],'orientation','vertical','AxisLocation','out')
%set(gcf,'renderer','painters')
text(0,1.75,'Mean fraction MZ hploss consumed','HorizontalAlignment','center')
text(-2.75,1.75,'A')

%2 - m over
subplot('Position',[0.475 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,HOmz*-1)
%colormap(cmBP2)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
%colorbar('Position',[0.825 0.695 0.025 0.275],'orientation','vertical','AxisLocation','out')
%set(gcf,'renderer','painters')
text(0,1.75,'Mean times MZ overconsumed','HorizontalAlignment','center')
text(-2.75,1.75,'B')

%3 - l frac
subplot('Position',[0.01 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,HFlz*-1)
%colormap(cmBP2)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
%colorbar('Position',[0.385 0.075 0.025 0.275],'orientation','vertical','AxisLocation','out')
%set(gcf,'renderer','painters')
text(0,1.75,'Mean fraction LZ hploss consumed','HorizontalAlignment','center')
text(-2.75,1.75,'C')

%4 - l over
subplot('Position',[0.475 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,HOlz*-1)
%colormap(cmBP2)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
%colorbar('Position',[0.825 0.075 0.025 0.275],'orientation','vertical','AxisLocation','out')
%set(gcf,'renderer','painters')
text(0,1.75,'Mean times MZ overconsumed','HorizontalAlignment','center')
text(-2.75,1.75,'D')
%print('-dpng',[ppath 'Hist_' harv '_global_zoop_overcon.png'])

%%
figure(2)
%1 - m frac
subplot('Position',[0.01 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,t1*-1)
%colormap(cmBP2)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
%colorbar('Position',[0.385 0.695 0.025 0.275],'orientation','vertical','AxisLocation','out')
%set(gcf,'renderer','painters')
text(0,1.75,'Mean fraction MZ hploss consumed','HorizontalAlignment','center')
text(-2.75,1.75,'A')

%2 - m over
subplot('Position',[0.475 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,t3*-1)
%colormap(cmBP2)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
%colorbar('Position',[0.825 0.695 0.025 0.275],'orientation','vertical','AxisLocation','out')
%set(gcf,'renderer','painters')
text(0,1.75,'Mean times MZ overconsumed','HorizontalAlignment','center')
text(-2.75,1.75,'B')

%3 - l frac
subplot('Position',[0.01 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,t4*-1)
%colormap(cmBP2)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
%colorbar('Position',[0.385 0.075 0.025 0.275],'orientation','vertical','AxisLocation','out')
%set(gcf,'renderer','painters')
text(0,1.75,'Mean fraction LZ hploss consumed','HorizontalAlignment','center')
text(-2.75,1.75,'C')

%4 - l over
subplot('Position',[0.475 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,t5*-1)
%colormap(cmBP2)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
%colorbar('Position',[0.825 0.075 0.025 0.275],'orientation','vertical','AxisLocation','out')
%set(gcf,'renderer','painters')
text(0,1.75,'Mean times MZ overconsumed','HorizontalAlignment','center')
text(-2.75,1.75,'D')
%print('-dpng',[ppath 'Hist_' harv '_global_zoop_overcon.png'])

%%
figure(10)
%1 - m frac
subplot('Position',[0.01 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([-2 2]);
colorbar
%set(gcf,'renderer','painters')
text(-2.75,1.75,'A')

%2 - m over
subplot('Position',[0.475 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllM))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([-2 2]);
colorbar
%set(gcf,'renderer','painters')
text(-2.75,1.75,'B')

%3 - l frac
subplot('Position',[0.01 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllL))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
text(-2.75,1.75,'C')

%4 - l over
subplot('Position',[0.475 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Pb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
text(-2.75,1.75,'D')
%print('-dpng',[ppath 'Hist_' harv '_global_zoop_overcon.png'])

%% individual plots
figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,HOmz*-1)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Plankton predictive skill')
print('-dpng',[ppath 'Plankton_skill_map.png'])
%print('-dpng',['Plankton_skill_map.png'])

figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,t5*-1)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('Fish predictive skill')
print('-dpng',[ppath 'Fish_skill_map.png'])
%print('-dpng',['Fish_skill_map.png'])

%%
figure(5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllM))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4.5 -0.5]);
colorbar('Ticks',[-4.5, -2.5, -0.5],'TickLabels',{'Low','Mid','High'})
set(gcf,'renderer','painters')
title('Predicted fish biomass')
print('-dpng',[ppath 'Fish_biomass_map.png'])

%% save
plank_skill = HOmz*-1;
fish_skill = t5*-1;
fish_mass = log10(AllM);

save([ppath 'skill_map_data.mat'],'geolat_t','geolon_t','plank_skill',...
    'fish_skill','fish_mass');

Mvec(:,1) = geolat_t(ID);
Mvec(:,2) = geolon_t(ID);
Mvec(:,3) = plank_skill(ID);
Mvec(:,4) = fish_skill(ID);
Mvec(:,5) = fish_mass(ID);
Mtab = array2table(Mvec,'VariableNames',{'Lat','Lon','PlankSkill','FishSkill',...
    'FishBiomass'});
writetable(Mtab,[ppath 'skill_map_data.csv'],'Delimiter',',')



