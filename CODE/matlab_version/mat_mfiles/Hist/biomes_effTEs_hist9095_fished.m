% Visualize output of FEISTY Historic globally
% 1990-1995, monthly means saved
% Transfer efficiency ("effective") 
% Use BE*det

clear all
close all

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
gpath='/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

% POEM
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
BE = 0.075;
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'TEeffDet_Historic9095_All_fish03_' cfile '.mat']);

cmYOR=cbrewer('seq','YlOrRd',50);
cmRP=cbrewer('seq','RdPu',50);
cmPR=cbrewer('seq','PuRd',50);

%% plot info

geolon_t=double(geolon_t);
geolat_t=double(geolat_t);
[ni,nj]=size(geolon_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

%% Effective TEs
% With BE*det instead of Bent

% Take means over biomes
load([gpath 'COBALT_biomes_ESM2M_9095.mat'])
tlme = biome_hist;

biome_te = NaN*ones(3,4);
for L=1:3
    lid = find(tlme==L);
    %TEeff
    biome_te(L,1) = nanmean(TEeffM(lid));
    biome_te(L,2) = nanmean(TEeff_L(lid));
    biome_te(L,3) = nanmean(TEeff_HTLd(lid));
    biome_te(L,4) = nanmean(TEeff_LTLd(lid));
    
end

biome_m = NaN*ones(ni,nj);
biome_l = biome_m;
biome_htlD = biome_m;
biome_ltlD = biome_m;
for L=1:3
    lid = find(tlme==L);

    biome_m(lid)      = biome_te(L,1);
    biome_l(lid)      = biome_te(L,2);
    biome_htlD(lid)   = biome_te(L,3);
    biome_ltlD(lid)   = biome_te(L,4);
end

%% Conversions to TE and save
TELTLd1 = real(biome_te(:,4).^(1/1.25));
TELTLd2 = real(biome_te(:,4).^(1/1.5));
TEM = real(biome_te(:,1).^(1/2));          %should this be 1/1?
TEL = real(biome_te(:,2).^(1/4));         %should this be 1/3?
TEHTLd = real(biome_te(:,3).^(1/3));   %should this be 1/2?

q(:,1) = TEM;
q(:,2) = TELTLd1;
q(:,3) = TELTLd2;
q(:,4) = TEHTLd;
q(:,5) = TEL;

Q = array2table(q,'VariableNames',{'TEM','TELTL1','TELTL2','TEHTL','TEATL'},...
    'RowNames',{'LC','ECCS','ECSS'});

B = array2table(biome_te,'VariableNames',{'TEM','TEATL','TEHTL','TELTL'},...
    'RowNames',{'LC','ECCS','ECSS'});

%% save
writetable(Q,[fpath 'TE_biomes_Historic9095_All_fish03_' cfile '.csv'],'Delimiter',',');
writetable(B,[fpath 'TEeff_biomes_Historic9095_All_fish03_' cfile '.csv'],'Delimiter',',');


%% Figures
% All 3 on subplots
%Detritus----------------------
figure(12)
subplot('Position',[0 0.53 0.5 0.5])
%LTL
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(TEeff_LTLd))
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 -0.6]);
colorbar('Position',[0.025 0.555 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('log_1_0 TEeff LTL')
text(-2.75,1.25,'A')

%HTL
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(TEeff_HTLd))
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3.5 -1]);
colorbar('Position',[0.525 0.555 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('log_1_0 TEeff HTL')
text(-2.75,1.25,'B')

%L
subplot('Position',[0.25 0.025 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(TEeff_L))
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5.5 -1.5]);
colorbar('Position',[0.275 0.05 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('log_1_0 TEeff L')
text(-2.75,1.25,'C')
%stamp([harv '_' cfile])
print('-dpng',[ppath 'Historic9095_' harv '_global_DeffTEs_subplot.png'])

%% All 3 converted on subplots
%Detritus----------------------
figure(14)
subplot('Position',[0 0.53 0.5 0.5])
%LTL
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(TEeff_LTLd))
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.03 0.15]);
colorbar('Position',[0.025 0.555 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('A. TE LTL')

%L
subplot('Position',[0.25 0.025 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TEL)
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.3]);
colorbar('Position',[0.275 0.05 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('C. TE L')

subplot('Position',[0.5 0.53 0.5 0.5])
%HTL
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,TEHTLd)
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.4]);
colorbar('Position',[0.525 0.555 0.45 0.05],'orientation','horizontal')                   
set(gcf,'renderer','painters')
title('B. TE HTL')
%stamp([harv '_' cfile])
print('-dpng',[ppath 'Historic9095_' harv '_global_DTEs_subplot.png'])


