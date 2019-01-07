% Calc LME biomass of POEM
% Climatology
% 150 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cdir 'temp_100_1deg_ESM26_5yr_clim_191_195.mat'])
load([cdir 'btm_temp_1deg_ESM26_5yr_clim_191_195.mat'])

ptemp_mean_clim=squeeze(nanmean(temp_100,1));
btemp_mean_clim=squeeze(nanmean(btm_temp,1));

%%
AREA_OCN = max(area,1);

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);
%load([dpath 'Means_Climatol_' harv '_' cfile '.mat']);

%% Plots in space
[ni,nj]=size(lon);

Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zb=NaN*ones(ni,nj);
Zd1=NaN*ones(ni,nj);
Zd2=NaN*ones(ni,nj);
Zm1=NaN*ones(ni,nj);
Zm2=NaN*ones(ni,nj);

Tsf=NaN*ones(ni,nj);
Tsp=NaN*ones(ni,nj);
Tsd=NaN*ones(ni,nj);
Tmf=NaN*ones(ni,nj);
Tmp=NaN*ones(ni,nj);
Tmd=NaN*ones(ni,nj);
Tlp=NaN*ones(ni,nj);
Tld=NaN*ones(ni,nj);

Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);

Zsf(ID)=sf_mean;
Zsp(ID)=sp_mean;
Zsd(ID)=sd_mean;
Zmf(ID)=mf_mean;
Zmp(ID)=mp_mean;
Zmd(ID)=md_mean;
Zlp(ID)=lp_mean;
Zld(ID)=ld_mean;
Zb(ID)=b_mean;
Zd1(ID)=all_median1;
Zd2(ID)=all_median2;
Zm1(ID)=all_mean1;
Zm2(ID)=all_mean2;

Tsf(ID)=sf_tot;
Tsp(ID)=sp_tot;
Tsd(ID)=sd_tot;
Tmf(ID)=mf_tot;
Tmp(ID)=mp_tot;
Tmd(ID)=md_tot;
Tlp(ID)=lp_tot;
Tld(ID)=ld_tot;

Cmf(ID)=mf_my;
Cmp(ID)=mp_my;
Cmd(ID)=md_my;
Clp(ID)=lp_my;
Cld(ID)=ld_my;

% g/m2/d --> total g
Amf_mcatch = Cmf .* AREA_OCN * 365; %mean fish catch per yr
Amp_mcatch = Cmp .* AREA_OCN * 365;
Amd_mcatch = Cmd .* AREA_OCN * 365;
Alp_mcatch = Clp .* AREA_OCN * 365;
Ald_mcatch = Cld .* AREA_OCN * 365;
% g/m2 --> total g
Asf_mean = Zsf .* AREA_OCN;
Asp_mean = Zsp .* AREA_OCN;
Asd_mean = Zsd .* AREA_OCN;
Amf_mean = Zmf .* AREA_OCN;
Amp_mean = Zmp .* AREA_OCN;
Amd_mean = Zmd .* AREA_OCN;
Alp_mean = Zlp .* AREA_OCN;
Ald_mean = Zld .* AREA_OCN;
Ab_mean  = Zb .* AREA_OCN;
all_men1 = Zm1 .* AREA_OCN;
all_men2 = Zm2 .* AREA_OCN;
all_med1 = Zd1 .* AREA_OCN;
all_med2 = Zd2 .* AREA_OCN;

Asf_tot = Tsf .* AREA_OCN;
Asp_tot = Tsp .* AREA_OCN;
Asd_tot = Tsd .* AREA_OCN;
Amf_tot = Tmf .* AREA_OCN;
Amp_tot = Tmp .* AREA_OCN;
Amd_tot = Tmd .* AREA_OCN;
Alp_tot = Tlp .* AREA_OCN;
Ald_tot = Tld .* AREA_OCN;

%% Total biomass (compare to J&C15)
tot_bio1 = nansum(all_med1(:)) * 1e-6
tot_bio2 = nansum(all_med2(:)) * 1e-6
tot_bio3 = nansum(all_men1(:)) * 1e-6
tot_bio4 = nansum(all_men2(:)) * 1e-6

%Too high by 1e2
tot_bio5 = nansum(Asf_tot(:) + Asp_tot(:) + Asd_tot(:) +...
    Amf_tot(:) + Amp_tot(:) + Amd_tot(:) +...
    Alp_tot(:) + Ald_tot(:)) * 1e-6;
tot_bio6 = nansum(Asf_mean(:) + Asp_mean(:) + Asd_mean(:) +...
    Amf_mean(:) + Amp_mean(:) + Amd_mean(:) +...
    Alp_mean(:) + Ald_mean(:)) * 1e-6;
tot_bioML = nansum(Amf_mean(:) + Amp_mean(:) + Amd_mean(:) +...
    Alp_mean(:) + Ald_mean(:)) * 1e-6;


%% Calc LMEs
tlme = lme_mask_onedeg;

lme_mcatch = NaN*ones(66,5);
lme_mbio = NaN*ones(66,9);
lme_sbio = NaN*ones(66,9);
lme_area = NaN*ones(66,1);
lme_med = NaN*ones(66,2);

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
    %median/mean of all
    lme_med(L,1) = nansum(all_med1(lid));
    lme_med(L,2) = nansum(all_men1(lid));
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

Tlme_med = sum(lme_med)* 1e-6;
frac_med_lme = Tlme_med(1)/tot_bio1
frac_men_lme = Tlme_med(2)/tot_bio3

%%
save([dpath 'LME_clim_fished_',harv,'_' cfile '.mat'],...
    'lme_mcatch','lme_mbio','lme_sbio','lme_area');

%% Figures

% clme_mf = NaN*ones(360,200);
% clme_mp = clme_mf;
% clme_md = clme_mf;
% clme_lp = clme_mf;
% clme_ld = clme_mf;
%
% lme_sf = NaN*ones(360,200);
% lme_sp = lme_sf;
% lme_sd = lme_sf;
% lme_mf = lme_sf;
% lme_mp = lme_sf;
% lme_md = lme_sf;
% lme_lp = lme_sf;
% lme_ld = lme_sf;
% lme_b = lme_sf;
%
% for L=1:66
%     lid = find(tlme==L);
%
%     clme_mf(lid) = lme_mcatch(L,1);
%     clme_mp(lid) = lme_mcatch(L,2);
%     clme_md(lid) = lme_mcatch(L,3);
%     clme_lp(lid) = lme_mcatch(L,4);
%     clme_ld(lid) = lme_mcatch(L,5);
%
%     lme_sf(lid) = lme_mbio(L,1);
%     lme_sp(lid) = lme_mbio(L,2);
%     lme_sd(lid) = lme_mbio(L,3);
%     lme_mf(lid) = lme_mbio(L,4);
%     lme_mp(lid) = lme_mbio(L,5);
%     lme_md(lid) = lme_mbio(L,6);
%     lme_lp(lid) = lme_mbio(L,7);
%     lme_ld(lid) = lme_mbio(L,8);
%     lme_b(lid) = lme_mbio(L,9);
% end
%
% clme_All = clme_mf+clme_mp+clme_md+clme_lp+clme_ld;
% clme_AllF = clme_mf;
% clme_AllP = clme_mp+clme_lp;
% clme_AllD = clme_md+clme_ld;
% clme_AllM = clme_mf+clme_mp+clme_md;
% clme_AllL = clme_lp+clme_ld;
%
% lme_All = lme_sf+lme_sp+lme_sd+lme_mf+lme_mp+lme_md+lme_lp+lme_ld;
% lme_AllF = lme_sf+lme_mf;
% lme_AllP = lme_sp+lme_mp+lme_lp;
% lme_AllD = lme_sd+lme_md+lme_ld;
%
% % plot info
% geolon_t = double(geolon_t);
% geolat_t = double(geolat_t);
% plotminlat=-90; %Set these bounds for your data
% plotmaxlat=90;
% plotminlon=-280;
% plotmaxlon=80;
% latlim=[plotminlat plotmaxlat];
% lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
% % ENTER -100 TO MAP ORIGIN LONG
%
% land=-999*ones(ni,nj);
% land(grid(:,1))=NaN*ones(size(mf_mean50));
%
% %% Catch
%
% % all
% figure(1)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_All*1e-6)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([4 7]);
% hcb = colorbar('h');
% ylim(hcb,[4 7])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Historic 1956-2005 LME mean log10 total annual catch (MT) All Fishes')
% stamp(cfile)
% print('-dpng',[ppath 'Clim_fished_',harv,'_LME_catch_All.png'])
%
% %% all F
% figure(2)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_AllF*1e-6)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([3.5 6.5]);
% hcb = colorbar('h');
% ylim(hcb,[3.5 6.5])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Historic 2051-2100 LME mean log10 total annual catch (MT) Forage Fishes')
% stamp(cfile)
% print('-dpng',[ppath 'Clim_fished_',harv,'_LME_catch_AllF.png'])
%
% % all P
% figure(3)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_AllP*1e-6)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([3.5 6.5]);
% hcb = colorbar('h');
% ylim(hcb,[3.5 6.5])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Historic 2051-2100 LME mean log10 total annual catch (MT) Pelagic Fishes')
% stamp(cfile)
% print('-dpng',[ppath 'Clim_fished_',harv,'_LME_catch_AllP.png'])
%
% % All D
% figure(4)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_AllD*1e-6)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([3.5 6.5]);
% hcb = colorbar('h');
% ylim(hcb,[3.5 6.5])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Historic 2051-2100 LME mean log10 total annual catch (MT) Demersal Fishes')
% stamp(cfile)
% print('-dpng',[ppath 'Clim_fished_',harv,'_LME_catch_AllD.png'])
%
% % all M
% figure(5)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_AllM*1e-6)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([3.5 6.5]);
% hcb = colorbar('h');
% ylim(hcb,[3.5 6.5])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Historic 2051-2100 LME mean log10 total annual catch (MT) Medium Fishes')
% stamp(cfile)
% print('-dpng',[ppath 'Clim_fished_',harv,'_LME_catch_AllM.png'])
%
% % all L
% figure(6)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_AllL*1e-6)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([3.5 6.5]);
% hcb = colorbar('h');
% ylim(hcb,[3.5 6.5])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Historic 2051-2100 LME mean log10 total annual catch (MT) Large Fishes')
% stamp(cfile)
% print('-dpng',[ppath 'Clim_fished_',harv,'_LME_catch_AllL.png'])
%



