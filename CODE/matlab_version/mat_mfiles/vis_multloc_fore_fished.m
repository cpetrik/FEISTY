% Visualize output of POEM
% Forecast time period (2006-2100) at all locations
% 95 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70_cmax-metab20_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00400';
harv = '03';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];

load([fpath 'Means_fore_fished',harv,'_' cfile '.mat']);

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);


%% Pick which time period mean
% 2051-2100
sp_smean=sp_mean50;
sf_smean=sf_mean50;
sd_smean=sd_mean50;
mp_smean=mp_mean50;
mf_smean=mf_mean50;
md_smean=md_mean50;
lp_smean=lp_mean50;
ld_smean=ld_mean50;
b_smean=b_mean50;

sp_sprod=sp_prod50;
sf_sprod=sf_prod50;
sd_sprod=sd_prod50;
mp_sprod=mp_prod50;
mf_sprod=mf_prod50;
md_sprod=md_prod50;
lp_sprod=lp_prod50;
ld_sprod=ld_prod50;

%% colors
cm9=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...  %b
    0 0 0];...      %black
    
cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm9);

%% Plots in time
y = 2005+(1/12):(1/12):2100;

% Piscivore
figure(1)
subplot(4,1,1)
plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Forecast Fished Pelagic Piscivores')
ylabel('log10 Biomass (g m^-^2)')
legend('Larvae','Juveniles','Adults')
legend('location','southeast')
stamp(cfile)

subplot(4,1,2)
plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Larvae')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,3)
plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Juveniles')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,4)
plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Adults')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_fished',harv,'_P_time.png'])

% Planktivore
sf_tmean=sf_tmean(1:length(y));
figure(2)
subplot(3,1,1)
plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Forecast Fished Forage Fishes')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
legend('Immature','Adults')
legend('location','southeast')
stamp(cfile)

subplot(3,1,2)
plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Immature')
ylabel('log10 Biomass (g m^-^2)')

subplot(3,1,3)
plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Adults')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_fished',harv,'_F_time.png'])

% Detritivore
figure(3)
subplot(4,1,1)
plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Forecast Fished Demersal Piscivores')
ylabel('log10 Biomass (g m^-^2)')
legend('Larvae','Juveniles','Adults')
legend('location','southeast')
stamp(cfile)

subplot(4,1,2)
plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Larvae')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,3)
plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Juveniles')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,4)
plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Adults')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_fished',harv,'_D_time.png'])

%% All size classes of all

figure(4)
plot(y,log10(sf_tmean),'Linewidth',1); hold on;
plot(y,log10(mf_tmean),'Linewidth',1); hold on;
plot(y,log10(sp_tmean),'Linewidth',1); hold on;
plot(y,log10(mp_tmean),'Linewidth',1); hold on;
plot(y,log10(lp_tmean),'Linewidth',1); hold on;
plot(y,log10(sd_tmean),'Linewidth',1); hold on;
plot(y,log10(md_tmean),'Linewidth',1); hold on;
plot(y,log10(ld_tmean),'Linewidth',1); hold on;
legend('SF','MF','SP','MP','LP','SD','MD','LD')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-5 2])
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
title('Forecast Fished')
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_all_sizes.png'])

figure(5)
F = sf_tmean+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;
D = sd_tmean+md_tmean+ld_tmean;

plot(y,log10(F),'r','Linewidth',2); hold on;
plot(y,log10(P),'b','Linewidth',2); hold on;
plot(y,log10(D),'k','Linewidth',2); hold on;
legend('F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-5 2])
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
title(['Forecast Fished'])
print('-dpng',[ppath 'Fore_fished',harv,'_all_types.png'])

%% FISHING All size classes of all

figure(6)
plot(y,log10(mf_tmy),'color',[0 0.7 0],'Linewidth',1); hold on;
plot(y,log10(mp_tmy),'color',[1 0 0],'Linewidth',1); hold on;
plot(y,log10(lp_tmy),'color',[0.5 0 0],'Linewidth',1); hold on;
plot(y,log10(md_tmy),'color',[0 0.5 0.75],'Linewidth',1); hold on;
plot(y,log10(ld_tmy),'color',[0 0 0.75],'Linewidth',1); hold on;
legend('MF','MP','LP','MD','LD')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-7 0])
xlabel('Time (mo)')
ylabel('log10 Catch (g m^-^2)')
title(['Fore fished',harv])
%stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv, '_catch_all_sizes.png'])

figure(7)
F = mf_tmy;
P = mp_tmy+lp_tmy;
D = md_tmy+ld_tmy;

plot(y,log10(F),'r','Linewidth',2); hold on;
plot(y,log10(P),'b','Linewidth',2); hold on;
plot(y,log10(D),'k','Linewidth',2); hold on;
legend('F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-7 0])
xlabel('Time (y)')
ylabel('log10 Catch (g m^-^2)')
title(['Fore fished',harv])
print('-dpng',[ppath 'Fore_fished',harv, '_catch_all_types.png'])

%% Recruitment
st=1:12:length(y);
en=12:12:length(y);
MPy = NaN*ones((length(time)/12),1);
SFy = MPy;
MDy = MPy;

for n=1:(length(time)/12)
    MPy(n) = nansum(mp_trec(st(n):en(n)));
    SFy(n) = nansum(sf_trec(st(n):en(n)));
    MDy(n) = nansum(md_trec(st(n):en(n)));
end

%%
figure(8)
subplot(3,1,2)
plot(2006:2100,log10(MPy),'b','Linewidth',2); hold on;
xlim([2006 2100])
ylabel('log10 annual Recruitment (g m^-^2)')
title('Pelagic piscivores')
stamp(cfile)

subplot(3,1,1)
plot(2006:2100,log10(SFy),'r','Linewidth',2); hold on;
xlim([2006 2100])
title('Forage fishes')

subplot(3,1,3)
plot(2006:2100,log10(MDy),'k','Linewidth',2); hold on;
xlim([2006 2100])
title('Demersal piscivores')
print('-dpng',[ppath 'Fore_fished',harv,'_recruitment.png'])

%% Time series
yr=y;

all_bio = sp_tmean+sf_tmean(1:length(time))+sd_tmean+mp_tmean+mf_tmean+md_tmean+lp_tmean+ld_tmean;

figure(70)
plot(yr,log10(all_bio),'k','LineWidth',2)
ylim([-2 2])
xlim([2005 2100])
xlabel('Year')
ylabel('All fish mean biomass (g/m^2)')
title('Forecast fished')
print('-dpng',[ppath 'Fore_fished',harv,'_ts_mbio.png'])

%% Plots in space
[ni,nj]=size(geolon_t);

Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zb=NaN*ones(ni,nj);

Zsf(grid(:,1))=sf_smean;
Zsp(grid(:,1))=sp_smean;
Zsd(grid(:,1))=sd_smean;
Zmf(grid(:,1))=mf_smean;
Zmp(grid(:,1))=mp_smean;
Zmd(grid(:,1))=md_smean;
Zlp(grid(:,1))=lp_smean;
Zld(grid(:,1))=ld_smean;
Zb(grid(:,1))=b_smean;

Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);

Cmf(grid(:,1))=mf_my50;
Cmp(grid(:,1))=mp_my50;
Cmd(grid(:,1))=md_my50;
Clp(grid(:,1))=lp_my50;
Cld(grid(:,1))=ld_my50;

%% bent
figure(9)
surf(geolon_t,geolat_t,log10(Zb)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean benthic biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2.5 0.5])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_BENT.png'])

%
mgZb = (Zb/9)*1e3;
figure(10)
surf(geolon_t,geolat_t,log10(mgZb)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean benthic biomass (mg C m^-^2)')
colormap('jet')
colorbar('h')
caxis([-0.8 2.3])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_BENT_mgC.png'])

% sp
figure(11)
surf(geolon_t,geolat_t,log10(Zsp)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Larval P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_SP.png'])

% sf
figure(12)
surf(geolon_t,geolat_t,log10(Zsf)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Larval F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_SF.png'])

% sd
figure(13)
surf(geolon_t,geolat_t,log10(Zsd)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Larval D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_SD.png'])

% mp
figure(14)
surf(geolon_t,geolat_t,log10(Zmp)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Juvenile P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_MP.png'])

% mf
figure(15)
surf(geolon_t,geolat_t,log10(Zmf)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Adult F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_MF.png'])

% md
figure(16)
surf(geolon_t,geolat_t,log10(Zmd)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Juvenile D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_MD.png'])

% lp
figure(17)
surf(geolon_t,geolat_t,log10(Zlp)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Adult P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_LP.png'])

% ld
figure(18)
surf(geolon_t,geolat_t,log10(Zld)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Adult D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_LD.png'])

%% Diff maps of all fish
All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllD = Zsd+Zmd+Zld;
AllS = Zsp+Zsf+Zsd;
AllM = Zmp+Zmf+Zmd;
AllL = Zlp+Zld;
FracPD = AllP ./ (AllP+AllD);
FracPF = AllP ./ (AllP+AllF);
FracPFvD = (AllP+AllF) ./ (AllP+AllF+AllD);

%% ALL
figure(21)
surf(geolon_t,geolat_t,log10(All)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean biomass All Fishes (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-1 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_All.png'])

% all F
figure(22)
surf(geolon_t,geolat_t,log10(AllF)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean biomass All F (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_AllF.png'])

% all D
figure(23)
surf(geolon_t,geolat_t,log10(AllD)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean biomass All D (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_AllD.png'])

% All P
figure(24)
surf(geolon_t,geolat_t,log10(AllP)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean biomass All P (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_AllP.png'])

% FracPD
figure(25)
surf(geolon_t,geolat_t,FracPD); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 P:D mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_FracPD.png'])

% FracPF
figure(26)
surf(geolon_t,geolat_t,FracPF); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 P:F mean biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_FracPF.png'])

%% Production
Psf=NaN*ones(ni,nj);
Psp=NaN*ones(ni,nj);
Psd=NaN*ones(ni,nj);
Pmf=NaN*ones(ni,nj);
Pmp=NaN*ones(ni,nj);
Pmd=NaN*ones(ni,nj);
Plp=NaN*ones(ni,nj);
Pld=NaN*ones(ni,nj);

Psf(grid(:,1))=sf_sprod;
Psp(grid(:,1))=sp_sprod;
Psd(grid(:,1))=sd_sprod;
Pmf(grid(:,1))=mf_sprod;
Pmp(grid(:,1))=mp_sprod;
Pmd(grid(:,1))=md_sprod;
Plp(grid(:,1))=lp_sprod;
Pld(grid(:,1))=ld_sprod;

%
Psp(Psp<=0)=NaN;
Psf(Psf<=0)=NaN;
Psd(Psd<=0)=NaN;
Pmp(Pmp<=0)=NaN;
Pmf(Pmf<=0)=NaN;
Pmd(Pmd<=0)=NaN;
Plp(Plp<=0)=NaN;
Pld(Pld<=0)=NaN;

%% sp
figure(31)
surf(geolon_t,geolat_t,log10(Psp)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Larval P production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_prod_SP.png'])

% sf
figure(32)
surf(geolon_t,geolat_t,log10(Psf)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Larval F production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_prod_SF.png'])

% sd
figure(33)
surf(geolon_t,geolat_t,log10(Psd)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Larval D production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_prod_SD.png'])

% mp
figure(34)
surf(geolon_t,geolat_t,log10(Pmp)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Juvenile P production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_prod_MP.png'])

% mf
figure(35)
surf(geolon_t,geolat_t,log10(Pmf)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Adult F production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_prod_MF.png'])

% md
figure(36)
surf(geolon_t,geolat_t,log10(Pmd)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Juvenile D production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_prod_MD.png'])

% lp
figure(37)
surf(geolon_t,geolat_t,log10(Plp)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Adult P production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_prod_LP.png'])

% ld
figure(38)
surf(geolon_t,geolat_t,log10(Pld)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean Adult D production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_prod_LD.png'])

% Diff maps of all fish
PAll = Psp+Psf+Psd+Pmp+Pmf+Pmd+Plp+Pld;
PAllF = Psf+Pmf;
PAllP = Psp+Pmp+Plp;
PAllD = Psd+Pmd+Pld;
PAllS = Psp+Psf+Psd;
PAllM = Pmp+Pmf+Pmd;
PAllL = Plp+Pld;

% ALL
figure(39)
surf(geolon_t,geolat_t,log10(PAll)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean production All Fishes (g m^-^2)')
colormap('jet')
colorbar('h')
%caxis([-1 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_prod_All.png'])

% all F
figure(40)
surf(geolon_t,geolat_t,log10(PAllF)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean production All F (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_prod_AllF.png'])

% all D
figure(41)
surf(geolon_t,geolat_t,log10(PAllD)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean production All D (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_prod_AllD.png'])

% All P
figure(42)
surf(geolon_t,geolat_t,log10(PAllP)); view(2); hold on;
shading flat
title('Forecast fished 2051-2100 log10 mean production All P (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv,'_global_prod_AllP.png'])

%% CATCH
% mp
figure(43)
surf(geolon_t,geolat_t,log10(Cmp)); view(2); hold on;
shading flat
title('log10 mean Juvenile P catch (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 -2])
%stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv, '_global_MP_catch.png'])

% mf
figure(44)
surf(geolon_t,geolat_t,log10(Cmf)); view(2); hold on;
shading flat
title('log10 mean Adult F catch (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 -2])
%stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv, '_global_MF_catch.png'])

% md
figure(45)
surf(geolon_t,geolat_t,log10(Cmd)); view(2); hold on;
shading flat
title('log10 mean Juvenile D catch (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 -2])
%stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv, '_global_MD_catch.png'])

% lp
figure(46)
surf(geolon_t,geolat_t,log10(Clp)); view(2); hold on;
shading flat
title('log10 mean Adult P catch (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 -2])
%stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv, '_global_LP_catch.png'])

% ld
figure(47)
surf(geolon_t,geolat_t,log10(Cld)); view(2); hold on;
shading flat
title('log10 mean Adult D catch (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 -2])
%stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv, '_global_LD_catch.png'])

%% Diff maps of all fish
CAll = Cmp+Cmf+Cmd+Clp+Cld;
CAllF = Cmf;
CAllP = Cmp+Clp;
CAllD = Cmd+Cld;
CAllM = Cmp+Cmf+Cmd;
CAllL = Clp+Cld;

% ALL
figure(48)
surf(geolon_t,geolat_t,log10(CAll)); view(2); hold on;
shading flat
title('log10 mean catch All Fishes (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 -2])
%stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv, '_global_All_catch.png'])

% all F
figure(49)
surf(geolon_t,geolat_t,log10(CAllF)); view(2); hold on;
shading flat
title('log10 mean catch All F (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 -2])
%stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv, '_global_AllF_catch.png'])

% all D
figure(50)
surf(geolon_t,geolat_t,log10(CAllD)); view(2); hold on;
shading flat
title('log10 mean catch All D (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 -2])
%stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv, '_global_AllD_catch.png'])

% All P
figure(51)
surf(geolon_t,geolat_t,log10(CAllP)); view(2); hold on;
shading flat
title('log10 mean catch All P (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 -2])
%stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv, '_global_AllP_catch.png'])

% all M
figure(52)
surf(geolon_t,geolat_t,log10(CAllM)); view(2); hold on;
shading flat
title('log10 mean catch All M (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 -2])
%stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv, '_global_AllM_catch.png'])

% All L
figure(53)
surf(geolon_t,geolat_t,log10(CAllL)); view(2); hold on;
shading flat
title('log10 mean catch All L (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 -2])
%stamp(cfile)
print('-dpng',[ppath 'Fore_fished',harv, '_global_AllL_catch.png'])
