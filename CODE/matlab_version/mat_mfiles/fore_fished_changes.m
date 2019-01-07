% Changes to global biomass 
% Changes to biomass in Pacific biomes

clear all
close all

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

cfile = 'Dc_enc70_cmax-metab20_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00400';
harv = '03';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([gpath 'hindcast_gridspec.mat'],'dat','geolat_t','geolon_t');
load([gpath 'lme_mask_esm2m.mat']);
load([gpath 'NPid_esm2m.mat']);
grid = csvread([gpath 'grid_csv.csv']);

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')
cmap1(1,:)=[1 1 1];
cmap1(2,:)=cmap_ppt(1,:);
cmap1(3,:)=cmap_ppt(3,:);
cmap1(4,:)=cmap_ppt(5,:);

cmap2(1,:)=cmap_ppt(1,:);
cmap2(2,:)=cmap_ppt(3,:);
cmap2(3,:)=cmap_ppt(5,:);

cmap3(1,:)=[0 0 0];
cmap3(2,:)=cmap_ppt(1,:);
cmap3(3,:)=cmap_ppt(3,:);
cmap3(4,:)=cmap_ppt(5,:);

cmap4(1,:)=cmap_ppt(3,:);
cmap4(2,:)=cmap_ppt(1,:);
cmap4(3,:)=cmap_ppt(5,:);

% figure
% y=ones(10,1);
% plot(1:10,y,'color',cmap_ppt(1,:)); hold on;    %blue
% plot(1:10,2*y,'color',cmap_ppt(2,:)); hold on;  %grey
% plot(1:10,3*y,'color',cmap_ppt(3,:)); hold on;  %orange
% plot(1:10,4*y,'color',cmap_ppt(4,:)); hold on;  %black
% plot(1:10,5*y,'color',cmap_ppt(5,:)); hold on;  %green
% plot(1:10,6*y,'color',cmap_ppt(6,:)); hold on;  %dark grey
% plot(1:10,7*y,'color',cmap_ppt(7,:)); hold on;  %light grey

% Biomes 1 = LC; 2 = ECCS; 3 = ECSS;

%% Hist prist
load([dpath 'LME_hist_pristine_' cfile '.mat'],'lme_mbio');
load([dpath 'Biomes_hist_pristine_' cfile '.mat'],'biome_mbio50');
load([dpath 'Means_hist_pristine_' cfile '.mat']);

histp_gm2 = lme_mbio;
histp_biome = biome_mbio50;

%Mean global biomass in g/m2
all_histp(:,1) = sf_mean50+sp_mean50+sd_mean50+mf_mean50+...
    mp_mean50+md_mean50+lp_mean50+ld_mean50;
all_histp(:,2) = sf_mean50+mf_mean50; %all f
all_histp(:,3) = sp_mean50+mp_mean50+lp_mean50; %all p
all_histp(:,4) = sd_mean50+md_mean50+ld_mean50; %all d
all_histp(:,5) = sf_mean50+sp_mean50+sd_mean50; %all s
all_histp(:,6) = mf_mean50+mp_mean50+md_mean50; %all m
all_histp(:,7) = lp_mean50+ld_mean50; %all l

%Time series
ts_histp(1,:) = sf_tmean;
ts_histp(2,:) = sp_tmean;
ts_histp(3,:) = sd_tmean;
ts_histp(4,:) = mf_tmean;
ts_histp(5,:) = mp_tmean;
ts_histp(6,:) = md_tmean;
ts_histp(7,:) = lp_tmean;
ts_histp(8,:) = ld_tmean;
ts_histp(9,:) = b_tmean;

%Mean total global biomass in g
tot_all_histp = all_histp .* repmat(grid(:,5),1,7);

clear lme_mbio biome_mbio50 sf_mean50 sp_mean50 sd_mean50
clear mf_mean50 mp_mean50 md_mean50 b_mean50 lp_mean50 ld_mean50
clear sf_mean sp_mean sd_mean mf_mean mp_mean md_mean b_mean lp_mean ld_mea
clear sf_tmean sp_tmean sd_tmean mf_tmean mp_tmean md_tmean b_tmean lp_tmean ld_tmean

%% Hist fished
load([dpath 'LME_hist_fished',harv,'_' cfile '.mat']);
load([dpath 'Biomes_hist_fished',harv,'_' cfile '.mat']);
load([dpath 'Means_hist_fished',harv,'_' cfile '.mat']);

histf_gm2 = lme_mbio;
histf_biome = biome_mbio50;

%Mean global biomass in g/m2
all_histf(:,1) = sf_mean50+sp_mean50+sd_mean50+mf_mean50+...
    mp_mean50+md_mean50+lp_mean50+ld_mean50;
all_histf(:,2) = sf_mean50+mf_mean50; %all f
all_histf(:,3) = sp_mean50+mp_mean50+lp_mean50; %all p
all_histf(:,4) = sd_mean50+md_mean50+ld_mean50; %all d
all_histf(:,5) = sf_mean50+sp_mean50+sd_mean50; %all s
all_histf(:,6) = mf_mean50+mp_mean50+md_mean50; %all m
all_histf(:,7) = lp_mean50+ld_mean50; %all l

%Time series
ts_histf(1,:) = sf_tmean;
ts_histf(2,:) = sp_tmean;
ts_histf(3,:) = sd_tmean;
ts_histf(4,:) = mf_tmean;
ts_histf(5,:) = mp_tmean;
ts_histf(6,:) = md_tmean;
ts_histf(7,:) = lp_tmean;
ts_histf(8,:) = ld_tmean;
ts_histf(9,:) = b_tmean;

ts_histf_catch(1,:) = mf_tmy;
ts_histf_catch(2,:) = mp_tmy;
ts_histf_catch(3,:) = md_tmy;
ts_histf_catch(4,:) = lp_tmy;
ts_histf_catch(5,:) = ld_tmy;

%Mean total global biomass in g
tot_all_histf = all_histf .* repmat(grid(:,5),1,7);

clear lme_mbio biome_mbio50 sf_mean50 sp_mean50 sd_mean50
clear mf_mean50 mp_mean50 md_mean50 b_mean50 lp_mean50 ld_mean50
clear sf_mean sp_mean sd_mean mf_mean mp_mean md_mean b_mean lp_mean ld_mea
clear sf_tmean sp_tmean sd_tmean mf_tmean mp_tmean md_tmean b_tmean lp_tmean ld_tmean
clear mf_tmy mp_tmy md_tmy lp_tmy ld_tmy

%% Fore pristine
load([dpath 'LME_fore_pristine_' cfile '.mat'],'lme_mbio');
load([dpath 'Biomes_fore_pristine_' cfile '.mat'],'biome_mbio50');
load([dpath 'Means_fore_pristine_' cfile '.mat']);

forep_gm2 = lme_mbio;
forep_biome = biome_mbio50;

%Mean global biomass in g/m2
all_forep(:,1) = sf_mean50+sp_mean50+sd_mean50+mf_mean50+...
    mp_mean50+md_mean50+lp_mean50+ld_mean50;
all_forep(:,2) = sf_mean50+mf_mean50; %all f
all_forep(:,3) = sp_mean50+mp_mean50+lp_mean50; %all p
all_forep(:,4) = sd_mean50+md_mean50+ld_mean50; %all d
all_forep(:,5) = sf_mean50+sp_mean50+sd_mean50; %all s
all_forep(:,6) = mf_mean50+mp_mean50+md_mean50; %all m
all_forep(:,7) = lp_mean50+ld_mean50; %all l

%Time series
ts_forep(1,:) = sf_tmean;
ts_forep(2,:) = sp_tmean;
ts_forep(3,:) = sd_tmean;
ts_forep(4,:) = mf_tmean;
ts_forep(5,:) = mp_tmean;
ts_forep(6,:) = md_tmean;
ts_forep(7,:) = lp_tmean;
ts_forep(8,:) = ld_tmean;
ts_forep(9,:) = b_tmean;

%Mean total global biomass in g
tot_all_forep = all_forep .* repmat(grid(:,5),1,7);

clear lme_mbio biome_mbio50 sf_mean50 sp_mean50 sd_mean50
clear mf_mean50 mp_mean50 md_mean50 b_mean50 lp_mean50 ld_mean50
clear sf_mean sp_mean sd_mean mf_mean mp_mean md_mean b_mean lp_mean ld_mea
clear sf_tmean sp_tmean sd_tmean mf_tmean mp_tmean md_tmean b_tmean lp_tmean ld_tmean

%% Fore fished
load([dpath 'LME_fore_fished',harv,'_' cfile '.mat'],'lme_mbio');
load([dpath 'Biomes_fore_fished',harv,'_' cfile '.mat'],'biome_mbio50');
load([dpath 'Means_fore_fished',harv,'_' cfile '.mat']);

foref_gm2 = lme_mbio;
foref_biome = biome_mbio50;

%Mean global biomass in g/m2
all_foref(:,1) = sf_mean50+sp_mean50+sd_mean50+mf_mean50+...
    mp_mean50+md_mean50+lp_mean50+ld_mean50;
all_foref(:,2) = sf_mean50+mf_mean50; %all f
all_foref(:,3) = sp_mean50+mp_mean50+lp_mean50; %all p
all_foref(:,4) = sd_mean50+md_mean50+ld_mean50; %all d
all_foref(:,5) = sf_mean50+sp_mean50+sd_mean50; %all s
all_foref(:,6) = mf_mean50+mp_mean50+md_mean50; %all m
all_foref(:,7) = lp_mean50+ld_mean50; %all l

%Time series
ts_foref(1,:) = sf_tmean;
ts_foref(2,:) = sp_tmean;
ts_foref(3,:) = sd_tmean;
ts_foref(4,:) = mf_tmean;
ts_foref(5,:) = mp_tmean;
ts_foref(6,:) = md_tmean;
ts_foref(7,:) = lp_tmean;
ts_foref(8,:) = ld_tmean;
ts_foref(9,:) = b_tmean;

ts_foref_catch(1,:) = mf_tmy;
ts_foref_catch(2,:) = mp_tmy;
ts_foref_catch(3,:) = md_tmy;
ts_foref_catch(4,:) = lp_tmy;
ts_foref_catch(5,:) = ld_tmy;

%Mean total global biomass in g
tot_all_foref = all_foref .* repmat(grid(:,5),1,7);

clear lme_mbio biome_mbio50 sf_mean50 sp_mean50 sd_mean50
clear mf_mean50 mp_mean50 md_mean50 b_mean50 lp_mean50 ld_mean50
clear sf_mean sp_mean sd_mean mf_mean mp_mean md_mean b_mean lp_mean ld_mea
clear sf_tmean sp_tmean sd_tmean mf_tmean mp_tmean md_tmean b_tmean lp_tmean ld_tmean
clear mf_tmy mp_tmy md_tmy lp_tmy ld_tmy

%% Global

Msum_histf = nansum(all_histf);
Msum_foref = nansum(all_foref);

Tsum_histf = nansum(tot_all_histf);
Tsum_foref = nansum(tot_all_foref);

%Future fishing vs. historic fishing
Mdiff_foref_histf = (Msum_foref - Msum_histf) ./ Msum_histf;
Tdiff_foref_histf = (Tsum_foref - Tsum_histf) ./ Tsum_histf;

%% Time series
ts_bio_histp = nansum(ts_histp(1:8,:));
Fts_bio_histp = (ts_histp(1,:)+ts_histp(4,:));
Pts_bio_histp = (ts_histp(2,:)+ts_histp(5,:)+ts_histp(7,:));
Dts_bio_histp = (ts_histp(3,:)+ts_histp(6,:)+ts_histp(8,:));

ts_bio_histf = nansum(ts_histf(1:8,:));
Fts_bio_histf = (ts_histf(1,:)+ts_histf(4,:));
Pts_bio_histf = (ts_histf(2,:)+ts_histf(5,:)+ts_histf(7,:));
Dts_bio_histf = (ts_histf(3,:)+ts_histf(6,:)+ts_histf(8,:));
ts_mcatch_histf = nansum(ts_histf_catch);
Fts_mcatch_histf = ts_histf_catch(1,:);
Pts_mcatch_histf = ts_histf_catch(2,:)+ts_histf_catch(4,:);
Dts_mcatch_histf = ts_histf_catch(3,:)+ts_histf_catch(5,:);

ts_bio_forep = nansum(ts_forep(1:8,:));
Fts_bio_forep = (ts_forep(1,:)+ts_forep(4,:));
Pts_bio_forep = (ts_forep(2,:)+ts_forep(5,:)+ts_forep(7,:));
Dts_bio_forep = (ts_forep(3,:)+ts_forep(6,:)+ts_forep(8,:));

ts_bio_foref = nansum(ts_foref(1:8,:));
Fts_bio_foref = (ts_foref(1,:)+ts_foref(4,:));
Pts_bio_foref = (ts_foref(2,:)+ts_foref(5,:)+ts_foref(7,:));
Dts_bio_foref = (ts_foref(3,:)+ts_foref(6,:)+ts_foref(8,:));
ts_mcatch_foref = nansum(ts_foref_catch);
Fts_mcatch_foref = ts_foref_catch(1,:);
Pts_mcatch_foref = ts_foref_catch(2,:)+ts_foref_catch(4,:);
Dts_mcatch_foref = ts_foref_catch(3,:)+ts_foref_catch(5,:);

%%
yrh=1860+(1/12):(1/12):2005;
yrf=2005+(1/12):(1/12):2100;

mFts_bio_histf = movmean(Fts_bio_histf,12);
mPts_bio_histf = movmean(Pts_bio_histf,12);
mDts_bio_histf = movmean(Dts_bio_histf,12);
mFts_bio_foref = movmean(Fts_bio_foref,12);
mPts_bio_foref = movmean(Pts_bio_foref,12);
mDts_bio_foref = movmean(Dts_bio_foref,12);

mFts_mcatch_histf = movmean(Fts_mcatch_histf,12);
mPts_mcatch_histf = movmean(Pts_mcatch_histf,12);
mDts_mcatch_histf = movmean(Dts_mcatch_histf,12);
mFts_mcatch_foref = movmean(Fts_mcatch_foref,12);
mPts_mcatch_foref = movmean(Pts_mcatch_foref,12);
mDts_mcatch_foref = movmean(Dts_mcatch_foref,12);

figure(1)
plot(yrh,mFts_bio_histf,'color',cmap_ppt(3,:),'LineWidth',2); hold on;
plot(yrh,mPts_bio_histf,'color',cmap_ppt(1,:),'LineWidth',2); hold on;
plot(yrh,mDts_bio_histf,'color',cmap_ppt(5,:),'LineWidth',2); hold on;
plot(yrf,mFts_bio_foref,'color',cmap_ppt(3,:),'LineWidth',2); hold on;
plot(yrf,mPts_bio_foref,'color',cmap_ppt(1,:),'LineWidth',2); hold on;
plot(yrf,mDts_bio_foref,'color',cmap_ppt(5,:),'LineWidth',2); hold on;
xlim([1860 2100])
xlabel('Year')
ylabel('Mean biomass (g/m^2)')
legend('F','P','D')
%legend('location','southwest')
title('Mean biomass')
print('-dpng',[ppath 'ts_mbio_type.png'])
%%
figure(2)
plot(yrh,mFts_mcatch_histf*365,'color',cmap_ppt(3,:),'LineWidth',2); hold on;
plot(yrh,mPts_mcatch_histf*365,'color',cmap_ppt(1,:),'LineWidth',2); hold on;
plot(yrh,mDts_mcatch_histf*365,'color',cmap_ppt(5,:),'LineWidth',2); hold on;
plot(yrf,mFts_mcatch_foref*365,'color',cmap_ppt(3,:),'LineWidth',2); hold on;
plot(yrf,mPts_mcatch_foref*365,'color',cmap_ppt(1,:),'LineWidth',2); hold on;
plot(yrf,mDts_mcatch_foref*365,'color',cmap_ppt(5,:),'LineWidth',2); hold on;
xlim([1860 2100])
%ylim([1500 6500])
xlabel('Year')
ylabel('Mean annual catch (g/m^2)')
legend('F','P','D')
%legend('location','southwest')
title('Mean catch')
print('-dpng',[ppath 'ts_mcatch_type.png'])



%% N Pac
[npid,ix,grid_np] = intersect(NPid,grid(:,1));

Msum_histf_np = nansum(all_histf(grid_np,:));
Msum_foref_np = nansum(all_foref(grid_np,:));

Tsum_histf_np = nansum(tot_all_histf(grid_np,:));
Tsum_foref_np = nansum(tot_all_foref(grid_np,:));

%Future fishing vs. historic fishing
Mdiff_foref_histf_np = (Msum_foref_np - Msum_histf_np) ./ Msum_histf_np;
Tdiff_foref_histf_np = (Tsum_foref_np - Tsum_histf_np) ./ Tsum_histf_np;

%% Maps
%plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;

x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);

%Global grid
All_histf_gg = griddata(grid(:,2),grid(:,3),all_histf(:,1),X,Y);
All_foref_gg = griddata(grid(:,2),grid(:,3),all_foref(:,1),X,Y);
F_histf_gg = griddata(grid(:,2),grid(:,3),all_histf(:,2),X,Y);
F_foref_gg = griddata(grid(:,2),grid(:,3),all_foref(:,2),X,Y);
P_histf_gg = griddata(grid(:,2),grid(:,3),all_histf(:,3),X,Y);
P_foref_gg = griddata(grid(:,2),grid(:,3),all_foref(:,3),X,Y);
D_histf_gg = griddata(grid(:,2),grid(:,3),all_histf(:,4),X,Y);
D_foref_gg = griddata(grid(:,2),grid(:,3),all_foref(:,4),X,Y);

%MOM grid
All_histf_mg = NaN*ones(size(geolon_t));
All_foref_mg = NaN*ones(size(geolon_t));
All_histf_mg(grid(:,1)) = all_histf(:,1);
All_foref_mg(grid(:,1)) = all_foref(:,1);

F_histf_mg = NaN*ones(size(geolon_t));
F_foref_mg = NaN*ones(size(geolon_t));
F_histf_mg(grid(:,1)) = all_histf(:,2);
F_foref_mg(grid(:,1)) = all_foref(:,2);

P_histf_mg = NaN*ones(size(geolon_t));
P_foref_mg = NaN*ones(size(geolon_t));
P_histf_mg(grid(:,1)) = all_histf(:,3);
P_foref_mg(grid(:,1)) = all_foref(:,3);

D_histf_mg = NaN*ones(size(geolon_t));
D_foref_mg = NaN*ones(size(geolon_t));
D_histf_mg(grid(:,1)) = all_histf(:,4);
D_foref_mg(grid(:,1)) = all_foref(:,4);

Diff_foref_histf_gg = (All_foref_gg - All_histf_gg) ./All_histf_gg;
Diff_foref_histf_mg = (All_foref_mg - All_histf_mg) ./All_histf_mg;
FDiff_foref_histf_gg = (F_foref_gg - F_histf_gg) ./F_histf_gg;
FDiff_foref_histf_mg = (F_foref_mg - F_histf_mg) ./F_histf_mg;
PDiff_foref_histf_gg = (P_foref_gg - P_histf_gg) ./P_histf_gg;
PDiff_foref_histf_mg = (P_foref_mg - P_histf_mg) ./P_histf_mg;
DDiff_foref_histf_gg = (D_foref_gg - D_histf_gg) ./D_histf_gg;
DDiff_foref_histf_mg = (D_foref_mg - D_histf_mg) ./D_histf_mg;

%% Global
%All
figure(3)
m_proj('miller','lat',82);
m_pcolor(X,Y,Diff_foref_histf_gg*100); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2051-2100 % difference from 1951-2000 mean biomass of all fishes')
colormap(cmap_color_rb)
colorbar('h')
caxis([-100 100])
stamp(cfile)
print('-dpng',[ppath 'Diff_global_all_foref_histf.png'])

%F
figure(4)
m_proj('miller','lat',82);
m_pcolor(X,Y,FDiff_foref_histf_gg*100); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2051-2100 % difference from 1951-2000 mean biomass of forage fishes')
colormap(cmap_color_rb)
colorbar('h')
caxis([-100 100])
stamp(cfile)
print('-dpng',[ppath 'Diff_global_F_foref_histf.png'])

%P
figure(5)
m_proj('miller','lat',82);
m_pcolor(X,Y,PDiff_foref_histf_gg*100); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2051-2100 % difference from 1951-2000 mean biomass of large pelagic fishes')
colormap(cmap_color_rb)
colorbar('h')
caxis([-100 100])
stamp(cfile)
print('-dpng',[ppath 'Diff_global_P_foref_histf.png'])

%D
figure(6)
m_proj('miller','lat',82);
m_pcolor(X,Y,DDiff_foref_histf_gg*100); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2051-2100 % difference from 1951-2000 mean biomass of demersal fishes')
colormap(cmap_color_rb)
colorbar('h')
caxis([-100 100])
stamp(cfile)
print('-dpng',[ppath 'Diff_global_D_foref_histf.png'])

%% N Pac
%All
figure(7)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,Diff_foref_histf_mg*100)
colormap(cmap_color_rb)              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
hcb = colorbar('h');
ylim(hcb,[-100 100])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference from 1951-2000 mean biomass of all fishes')
stamp(cfile)
print('-dpng',[ppath 'Diff_NPac_all_foref_histf.png'])

%F
figure(8)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,FDiff_foref_histf_mg*100)
colormap(cmap_color_rb)              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
hcb = colorbar('h');
ylim(hcb,[-100 100])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference from 1951-2000 mean biomass of forage fishes')
stamp(cfile)
print('-dpng',[ppath 'Diff_NPac_F_foref_histf.png'])

%P
figure(9)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,PDiff_foref_histf_mg*100)
colormap(cmap_color_rb)              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
hcb = colorbar('h');
ylim(hcb,[-100 100])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference from 1951-2000 mean biomass of large pelagic fishes')
stamp(cfile)
print('-dpng',[ppath 'Diff_NPac_P_foref_histf.png'])

%D
figure(10)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,DDiff_foref_histf_mg*100)
colormap(cmap_color_rb)              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
hcb = colorbar('h');
ylim(hcb,[-100 100])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference from 1951-2000 mean biomass of demersal fishes')
stamp(cfile)
print('-dpng',[ppath 'Diff_NPac_D_foref_histf.png'])




