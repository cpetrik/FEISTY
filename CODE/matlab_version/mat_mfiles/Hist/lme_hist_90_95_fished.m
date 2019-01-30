% Calc LME biomass of FEISTY
% Hindcast time period 1990-1995 (for comparing with Climatology)
% 5 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t','AREA_OCN');
grid = csvread([cpath 'grid_csv.csv']);
load([cpath 'lme_mask_esm2m.mat']);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

%% FEISTY
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'Means_Historic_',harv,'_' cfile '.mat']);

%% Plots in space
[ni,nj]=size(geolon_t);

Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);

Zmf(grid(:,1))=mf_my5;
Zmp(grid(:,1))=mp_my5;
Zmd(grid(:,1))=md_my5;
Zlp(grid(:,1))=lp_my5;
Zld(grid(:,1))=ld_my5;

Bsf = Zmf;
Bsp = Zmf;
Bsd = Zmf;
Bmf = Zmf;
Bmp = Zmf;
Bmd = Zmf;
Blp = Zmf;
Bld = Zmf;
Bb = Zmf;

Bsf(grid(:,1))=sf_mean5;
Bsp(grid(:,1))=sp_mean5;
Bsd(grid(:,1))=sd_mean5;
Bmf(grid(:,1))=mf_mean5;
Bmp(grid(:,1))=mp_mean5;
Bmd(grid(:,1))=md_mean5;
Blp(grid(:,1))=lp_mean5;
Bld(grid(:,1))=ld_mean5;
Bb(grid(:,1))=b_mean5;

% g/m2 --> total g
Amf_mcatch = Zmf .* AREA_OCN * 365; %mean fish catch per yr
Amp_mcatch = Zmp .* AREA_OCN * 365;
Amd_mcatch = Zmd .* AREA_OCN * 365;
Alp_mcatch = Zlp .* AREA_OCN * 365;
Ald_mcatch= Zld .* AREA_OCN * 365;

Asf_mean = Bsf .* AREA_OCN;
Asp_mean = Bsp .* AREA_OCN;
Asd_mean = Bsd .* AREA_OCN;
Amf_mean = Bmf .* AREA_OCN;
Amp_mean = Bmp .* AREA_OCN;
Amd_mean = Bmd .* AREA_OCN;
Alp_mean = Blp .* AREA_OCN;
Ald_mean = Bld .* AREA_OCN;
Ab_mean = Bb .* AREA_OCN;

%% Calc LMEs
lat = geolat_t;
lon = geolon_t;

tlme = lme_mask_esm2m';

lme_mcatch = NaN*ones(66,5);
lme_mbio = NaN*ones(66,9);
lme_sbio = NaN*ones(66,9);
lme_area = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    %total catch g
    lme_mcatch(L,1) = nansum(Amf_mcatch(lid));
    lme_mcatch(L,2) = nansum(Amp_mcatch(lid));
    lme_mcatch(L,3) = nansum(Amd_mcatch(lid));
    lme_mcatch(L,4) = nansum(Alp_mcatch(lid));
    lme_mcatch(L,5) = nansum(Ald_mcatch(lid));
    %mean biomass
    lme_mbio(L,1) = nanmean(Asf_mean(lid));
    lme_mbio(L,2) = nanmean(Asp_mean(lid));
    lme_mbio(L,3) = nanmean(Asd_mean(lid));
    lme_mbio(L,4) = nanmean(Amf_mean(lid));
    lme_mbio(L,5) = nanmean(Amp_mean(lid));
    lme_mbio(L,6) = nanmean(Amd_mean(lid));
    lme_mbio(L,7) = nanmean(Alp_mean(lid));
    lme_mbio(L,8) = nanmean(Ald_mean(lid));
    lme_mbio(L,9) = nanmean(Ab_mean(lid));
    %total biomass
    lme_sbio(L,1) = nansum(Asf_mean(lid));
    lme_sbio(L,2) = nansum(Asp_mean(lid));
    lme_sbio(L,3) = nansum(Asd_mean(lid));
    lme_sbio(L,4) = nansum(Amf_mean(lid));
    lme_sbio(L,5) = nansum(Amp_mean(lid));
    lme_sbio(L,6) = nansum(Amd_mean(lid));
    lme_sbio(L,7) = nansum(Alp_mean(lid));
    lme_sbio(L,8) = nansum(Ald_mean(lid));
    lme_sbio(L,9) = nansum(Ab_mean(lid));
    %total area of LME
    lme_area(L,1) = nansum(AREA_OCN(lid));
end

%%
save([dpath 'LME_hist_90-95_fished_',harv,'_' cfile '.mat'],...
    'lme_mcatch','lme_mbio','lme_sbio','lme_area');

%% Figures

clme_mf = NaN*ones(360,200);
clme_mp = clme_mf;
clme_md = clme_mf;
clme_lp = clme_mf;
clme_ld = clme_mf;

lme_sf = NaN*ones(360,200);
lme_sp = lme_sf;
lme_sd = lme_sf;
lme_mf = lme_sf;
lme_mp = lme_sf;
lme_md = lme_sf;
lme_lp = lme_sf;
lme_ld = lme_sf;
lme_b = lme_sf;

for L=1:66
    lid = find(tlme==L);
    
    clme_mf(lid) = lme_mcatch(L,1);
    clme_mp(lid) = lme_mcatch(L,2);
    clme_md(lid) = lme_mcatch(L,3);
    clme_lp(lid) = lme_mcatch(L,4);
    clme_ld(lid) = lme_mcatch(L,5);
    
    lme_sf(lid) = lme_mbio(L,1);
    lme_sp(lid) = lme_mbio(L,2);
    lme_sd(lid) = lme_mbio(L,3);
    lme_mf(lid) = lme_mbio(L,4);
    lme_mp(lid) = lme_mbio(L,5);
    lme_md(lid) = lme_mbio(L,6);
    lme_lp(lid) = lme_mbio(L,7);
    lme_ld(lid) = lme_mbio(L,8);
    lme_b(lid) = lme_mbio(L,9);
end

clme_All = clme_mf+clme_mp+clme_md+clme_lp+clme_ld;
clme_AllF = clme_mf;
clme_AllP = clme_mp+clme_lp;
clme_AllD = clme_md+clme_ld;
clme_AllM = clme_mf+clme_mp+clme_md;
clme_AllL = clme_lp+clme_ld;

lme_All = lme_sf+lme_sp+lme_sd+lme_mf+lme_mp+lme_md+lme_lp+lme_ld;
lme_AllF = lme_sf+lme_mf;
lme_AllP = lme_sp+lme_mp+lme_lp;
lme_AllD = lme_sd+lme_md+lme_ld;

% plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
% ENTER -100 TO MAP ORIGIN LONG

land=-999*ones(ni,nj);
land(grid(:,1))=NaN*ones(size(mf_mean5));

%% Catch

% all
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(clme_All*1e-6)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3.5 7.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Historic 1990-1994 LME mean log10 total annual catch (MT) All Fishes')
stamp(cfile)
%print('-dpng',[ppath 'Hist_fished',harv,'_LME_catch_All.png'])

%% all F
figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(clme_AllF*1e-6)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3.5 7.5]);
hcb = colorbar('h');
ylim(hcb,[3.5 7.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic 1990-1994 LME mean log10 total annual catch (MT) Forage Fishes')
stamp(cfile)
%print('-dpng',[ppath 'Hist_fished',harv,'_LME_catch_AllF.png'])

% all P
figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(clme_AllP*1e-6)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3.5 7.5]);
hcb = colorbar('h');
ylim(hcb,[3.5 7.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic 1990-1994 LME mean log10 total annual catch (MT) Pelagic Fishes')
stamp(cfile)
%print('-dpng',[ppath 'Hist_fished',harv,'_LME_catch_AllP.png'])

% All D
figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(clme_AllD*1e-6)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3.5 7.5]);
hcb = colorbar('h');
ylim(hcb,[3.5 7.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic 1990-1994 LME mean log10 total annual catch (MT) Demersal Fishes')
stamp(cfile)
%print('-dpng',[ppath 'Hist_fished',harv,'_LME_catch_AllD.png'])

% all M
figure(5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(clme_AllM*1e-6)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3.5 7.5]);
hcb = colorbar('h');
ylim(hcb,[3.5 7.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic 1990-1994 LME mean log10 total annual catch (MT) Medium Fishes')
stamp(cfile)
%print('-dpng',[ppath 'Hist_fished',harv,'_LME_catch_AllM.png'])

% all L
figure(6)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,real(log10(clme_AllL*1e-6)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3.5 7.5]);
hcb = colorbar('h');
ylim(hcb,[3.5 7.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Historic 1990-1994 LME mean log10 total annual catch (MT) Large Fishes')
stamp(cfile)
%print('-dpng',[ppath 'Hist_fished',harv,'_LME_catch_AllL.png'])

