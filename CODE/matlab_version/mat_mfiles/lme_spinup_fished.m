% Calc LME biomass of POEM
% Spinup globally
% 50 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

cfile = 'Diff_Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE10_CC100_RE00100';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'Means_spinup_fished_' cfile '.mat']);

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t','AREA_OCN');
grid = csvread([cpath 'grid_csv.csv']);
load([cpath 'lme_mask_esm2m.mat']);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

%% Plots in space
[ni,nj]=size(geolon_t);

Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);

Zmf(grid(:,1))=mf_my;
Zmp(grid(:,1))=mp_my;
Zmd(grid(:,1))=md_my;
Zlp(grid(:,1))=lp_my;
Zld(grid(:,1))=ld_my;

% g/m2 --> total g
Amf_mcatch = Zmf .* AREA_OCN * 365; %mean fish catch per yr
Amp_mcatch = Zmp .* AREA_OCN * 365; %mean fish catch per yr
Amd_mcatch = Zmd .* AREA_OCN * 365; %mean fish catch per yr
Alp_mcatch = Zlp .* AREA_OCN * 365;
Ald_mcatch= Zld .* AREA_OCN * 365;

%% Calc LMEs
lat = geolat_t;
lon = geolon_t;

tlme = lme_mask_esm2m';

lme_mcatch = NaN*ones(66,5);
lme_area = NaN*ones(66,1);
    
for L=1:66
    lid = find(tlme==L);
    %total catch g
    lme_mcatch(L,1) = nansum(Amf_mcatch(lid));
    lme_mcatch(L,2) = nansum(Amp_mcatch(lid));
    lme_mcatch(L,3) = nansum(Amd_mcatch(lid));
    lme_mcatch(L,4) = nansum(Alp_mcatch(lid));
    lme_mcatch(L,5) = nansum(Ald_mcatch(lid));
    lme_area(L,1) = nansum(AREA_OCN(lid));
end

%%
save([dpath 'LME_spinup_fished_' cfile '.mat'],'lme_area',...
    'lme_mcatch');


%% Figures

clme_mf = NaN*ones(360,200);
clme_mp = NaN*ones(360,200);
clme_md = NaN*ones(360,200);
clme_lp = clme_mf;
clme_ld = clme_mf;

for L=1:66
    lid = find(tlme==L);
    
    clme_mf(lid) = lme_mcatch(L);
    clme_mp(lid) = mp_lme_mcatch(L);
    clme_md(lid) = md_lme_mcatch(L);
    clme_lp(lid) = lp_lme_mcatch(L);
    clme_ld(lid) = ld_lme_mcatch(L);
end

clme_All = clme_mf+clme_mp+clme_md+clme_lp+clme_ld;
clme_AllF = clme_mf;
clme_AllP = clme_mp+clme_lp;
clme_AllD = clme_md+clme_ld;
clme_AllM = clme_mf+clme_mp+clme_md;
clme_AllL = clme_lp+clme_ld;

% plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
% ENTER -100 TO MAP ORIGIN LONG

land=-999*ones(ni,nj);
land(grid(:,1))=NaN*ones(size(mf_mean));
%% Catch

% all
figure(1)
surf(geolon_t,geolat_t,log10(clme_All*1e-6)); view(2); hold on;
shading flat
title('Spinup LME mean log10 total annual catch (MT) All Fishes')
colormap('jet')
colorbar('h')
caxis([4 7])
surf(geolon_t,geolat_t,land); view(2); hold on;
shading flat
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_LME_catch_All1.png'])

figure(51)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(clme_All*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([4 7]);
hcb = colorbar('h');
ylim(hcb,[4 7])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Spinup LME mean log10 total annual catch (MT) All Fishes')
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_LME_catch_All2.png'])

%% all F
figure(2)
surf(geolon_t,geolat_t,log10(clme_AllF*1e-6)); view(2); hold on;
shading flat
title('Spinup LME mean log10 total annual catch (MT) Forage Fishes')
colormap('jet')
colorbar('h')
caxis([3.5 6.5])
surf(geolon_t,geolat_t,land); view(2); hold on;
shading flat
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_LME_catch_AllF1.png'])

figure(52)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(clme_AllF*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3.5 6.5]);
hcb = colorbar('h');
ylim(hcb,[3.5 6.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Spinup LME mean log10 total annual catch (MT) Forage Fishes')
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_LME_catch_AllF2.png'])

%% all P
figure(3)
surf(geolon_t,geolat_t,log10(clme_AllP*1e-6)); view(2); hold on;
shading flat
title('Spinup LME mean log10 total annual catch (MT) Pelagic Fishes')
colormap('jet')
colorbar('h')
caxis([3.5 6.5])
surf(geolon_t,geolat_t,land); view(2); hold on;
shading flat
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_LME_catch_AllP1.png'])

figure(53)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(clme_AllP*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3.5 6.5]);
hcb = colorbar('h');
ylim(hcb,[3.5 6.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Spinup LME mean log10 total annual catch (MT) Pelagic Fishes')
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_LME_catch_AllP2.png'])

%% All D
figure(4)
surf(geolon_t,geolat_t,log10(clme_AllD*1e-6)); view(2); hold on;
shading flat
title('Spinup LME mean log10 total annual catch (MT) Demersal Fishes')
colormap('jet')
colorbar('h')
caxis([3.5 6.5])
surf(geolon_t,geolat_t,land); view(2); hold on;
shading flat
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_LME_catch_AllD1.png'])

figure(54)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(clme_AllD*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3.5 6.5]);
hcb = colorbar('h');
ylim(hcb,[3.5 6.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Spinup LME mean log10 total annual catch (MT) Demersal Fishes')
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_LME_catch_AllD2.png'])

% all M
figure(55)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(clme_AllM*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3.5 6.5]);
hcb = colorbar('h');
ylim(hcb,[3.5 6.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Spinup LME mean log10 total annual catch (MT) Medium Fishes')
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_LME_catch_AllM.png'])

%% all L
figure(6)
surf(geolon_t,geolat_t,log10(clme_AllL*1e-6)); view(2); hold on;
shading flat
title('Spinup LME mean log10 total annual catch (MT) Large Fishes')
colormap('jet')
colorbar('h')
caxis([3.5 6.5])
surf(geolon_t,geolat_t,land); view(2); hold on;
shading flat
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_LME_catch_AllL1.png'])

figure(56)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(clme_AllL*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3.5 6.5]);
hcb = colorbar('h');
ylim(hcb,[3.5 6.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Spinup LME mean log10 total annual catch (MT) Large Fishes')
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_LME_catch_AllL2.png'])

