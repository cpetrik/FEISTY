% Changes to global biomass 
% Changes to biomass in Pacific biomes

clear all
close all

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
dp = '/Volumes/GFDL/NC/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([gpath 'hindcast_gridspec.mat'],'dat','geolat_t','geolon_t');
load([gpath 'lme_mask_esm2m.mat']);
load([gpath 'NPid_esm2m.mat']);
load([cpath 'COBALT_biomes_last50yr.mat']);
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

%% Preindust
load([dpath 'LME_preindust_' cfile '.mat']);
load([dpath 'Biomes_preindust_' cfile '.mat'],'biome_mbio');
%load([dpath 'Preindust_1800-1850_means.mat']); %fcrit30
load([dpath 'Means_preindust_' cfile '.mat']); %fcrit40

pre_g = lme_bio;
pre_gm2 = lme_mbio;

pre_biome = biome_mbio;

%Mean global biomass in g/m2
all_pre(:,1) = sf_smean+sp_smean+sd_smean+mf_smean+mp_smean+md_smean+...
    lp_smean+ld_smean;
all_pre(:,2) = sf_smean+mf_smean; %all f
all_pre(:,3) = sp_smean+mp_smean+lp_smean; %all p
all_pre(:,4) = sd_smean+md_smean+ld_smean; %all d
all_pre(:,5) = sf_smean+sp_smean+sd_smean; %all s
all_pre(:,6) = mf_smean+mp_smean+md_smean; %all m
all_pre(:,7) = lp_smean+ld_smean; %all l

%Mean total global biomass in g
tot_all_pre = all_pre .* repmat(grid(:,5),1,7);

clear lme_bio lme_mbio biome_mbio sf_smean sp_smean sd_smean
clear mf_smean mp_smean md_smean b_smean lp_smean ld_smean

%% Hist prist
load([dpath 'LME_hist_pristine_' cfile '.mat'],'lme_bio00','lme_mbio00');
load([dpath 'Biomes_hist_pristine_' cfile '.mat'],'biome_mbio50');
load([dpath 'Means_hist_pristine_' cfile '.mat'],'sf_mean5000','sp_mean5000',...
    'sd_mean5000','mf_mean5000','mp_mean5000','md_mean5000','b_mean5000',...
    'lp_mean5000','ld_mean5000');

histp_g = lme_bio00;
histp_gm2 = lme_mbio00;

histp_biome = biome_mbio50;

%Mean global biomass in g/m2
all_histp(:,1) = sf_mean5000+sp_mean5000+sd_mean5000+mf_mean5000+...
    mp_mean5000+md_mean5000+lp_mean5000+ld_mean5000;
all_histp(:,2) = sf_mean5000+mf_mean5000; %all f
all_histp(:,3) = sp_mean5000+mp_mean5000+lp_mean5000; %all p
all_histp(:,4) = sd_mean5000+md_mean5000+ld_mean5000; %all d
all_histp(:,5) = sf_mean5000+sp_mean5000+sd_mean5000; %all s
all_histp(:,6) = mf_mean5000+mp_mean5000+md_mean5000; %all m
all_histp(:,7) = lp_mean5000+ld_mean5000; %all l

%Mean total global biomass in g
tot_all_histp = all_histp .* repmat(grid(:,5),1,7);

clear lme_bio00 lme_mbio00 biome_mbio50 sf_mean5000 sp_mean5000 sd_mean5000
clear mf_mean5000 mp_mean5000 md_mean5000 b_mean5000 lp_mean5000 ld_mean5000

%% Hist fished
load([dpath 'LME_hist_fished_' cfile '.mat'],'lme_bio00','lme_mbio00');
load([dpath 'Biomes_hist_fished_' cfile '.mat'],'biome_mbio50');
load([dpath 'Means_hist_fished_' cfile '.mat'],'sf_mean5000','sp_mean5000',...
    'sd_mean5000','mf_mean5000','mp_mean5000','md_mean5000','b_mean5000',...
    'lp_mean5000','ld_mean5000');

histf_g = lme_bio00;
histf_gm2 = lme_mbio00;

histf_biome = biome_mbio50;

%Mean global biomass in g/m2
all_histf(:,1) = sf_mean5000+sp_mean5000+sd_mean5000+mf_mean5000+...
    mp_mean5000+md_mean5000+lp_mean5000+ld_mean5000;
all_histf(:,2) = sf_mean5000+mf_mean5000; %all f
all_histf(:,3) = sp_mean5000+mp_mean5000+lp_mean5000; %all p
all_histf(:,4) = sd_mean5000+md_mean5000+ld_mean5000; %all d
all_histf(:,5) = sf_mean5000+sp_mean5000+sd_mean5000; %all s
all_histf(:,6) = mf_mean5000+mp_mean5000+md_mean5000; %all m
all_histf(:,7) = lp_mean5000+ld_mean5000; %all l

%Mean total global biomass in g
tot_all_histf = all_histf .* repmat(grid(:,5),1,7);

clear lme_bio00 lme_mbio00 biome_mbio50 sf_mean5000 sp_mean5000 sd_mean5000
clear mf_mean5000 mp_mean5000 md_mean5000 b_mean5000 lp_mean5000 ld_mean5000

%% Fore pristine
load([dpath 'LME_fore_pristine_' cfile '.mat'],'lme_bio00','lme_mbio00');
load([dpath 'Biomes_fore_pristine_' cfile '.mat'],'biome_mbio50');
load([dpath 'Means_fore_pristine_' cfile '.mat'],'sf_mean5000','sp_mean5000',...
    'sd_mean5000','mf_mean5000','mp_mean5000','md_mean5000','b_mean5000',...
    'lp_mean5000','ld_mean5000');

forep_g = lme_bio00;
forep_gm2 = lme_mbio00;

forep_biome = biome_mbio50;

%Mean global biomass in g/m2
all_forep(:,1) = sf_mean5000+sp_mean5000+sd_mean5000+mf_mean5000+...
    mp_mean5000+md_mean5000+lp_mean5000+ld_mean5000;
all_forep(:,2) = sf_mean5000+mf_mean5000; %all f
all_forep(:,3) = sp_mean5000+mp_mean5000+lp_mean5000; %all p
all_forep(:,4) = sd_mean5000+md_mean5000+ld_mean5000; %all d
all_forep(:,5) = sf_mean5000+sp_mean5000+sd_mean5000; %all s
all_forep(:,6) = mf_mean5000+mp_mean5000+md_mean5000; %all m
all_forep(:,7) = lp_mean5000+ld_mean5000; %all l

%Mean total global biomass in g
tot_all_forep = all_forep .* repmat(grid(:,5),1,7);

clear lme_bio00 lme_mbio00 biome_mbio50 sf_mean5000 sp_mean5000 sd_mean5000
clear mf_mean5000 mp_mean5000 md_mean5000 b_mean5000 lp_mean5000 ld_mean5000

%% Fore fished
load([dpath 'LME_fore_fished_' cfile '.mat'],'lme_bio00','lme_mbio00');
load([dpath 'Biomes_fore_fished_' cfile '.mat'],'biome_mbio50');
load([dpath 'Means_fore_fished_' cfile '.mat'],'sf_mean5000','sp_mean5000',...
    'sd_mean5000','mf_mean5000','mp_mean5000','md_mean5000','b_mean5000',...
    'lp_mean5000','ld_mean5000');

foref_g = lme_bio00;
foref_gm2 = lme_mbio00;

foref_biome = biome_mbio50;

%Mean global biomass in g/m2
all_foref(:,1) = sf_mean5000+sp_mean5000+sd_mean5000+mf_mean5000+...
    mp_mean5000+md_mean5000+lp_mean5000+ld_mean5000;
all_foref(:,2) = sf_mean5000+mf_mean5000; %all f
all_foref(:,3) = sp_mean5000+mp_mean5000+lp_mean5000; %all p
all_foref(:,4) = sd_mean5000+md_mean5000+ld_mean5000; %all d
all_foref(:,5) = sf_mean5000+sp_mean5000+sd_mean5000; %all s
all_foref(:,6) = mf_mean5000+mp_mean5000+md_mean5000; %all m
all_foref(:,7) = lp_mean5000+ld_mean5000; %all l

%Mean total global biomass in g
tot_all_foref = all_foref .* repmat(grid(:,5),1,7);

clear lme_bio00 lme_mbio00 biome_mbio50 sf_mean5000 sp_mean5000 sd_mean5000
clear mf_mean5000 mp_mean5000 md_mean5000 b_mean5000 lp_mean5000 ld_mean5000

%% Global

Msum_pre = nansum(all_pre);
Msum_histp = nansum(all_histp);
Msum_histf = nansum(all_histf);
Msum_forep = nansum(all_forep);
Msum_foref = nansum(all_foref);

Tsum_pre = nansum(tot_all_pre);
Tsum_histp = nansum(tot_all_histp);
Tsum_histf = nansum(tot_all_histf);
Tsum_forep = nansum(tot_all_forep);
Tsum_foref = nansum(tot_all_foref);

%Indust vs. Preindust
Mdiff_histp_pre = (Msum_histp - Msum_pre) ./ Msum_pre;
Tdiff_histp_pre = (Tsum_histp - Tsum_pre) ./ Tsum_pre;

%Hist fishing vs. PreIndust
Mdiff_histf_pre = (Msum_histf - Msum_pre) ./ Msum_pre;
Tdiff_histf_pre = (Tsum_histf - Tsum_pre) ./ Tsum_pre;

%Hist fishing vs. no fishing
Mdiff_histf_histp = (Msum_histf - Msum_histp) ./ Msum_histp;
Tdiff_histf_histp = (Tsum_histf - Tsum_histp) ./ Tsum_histp;

%Future fishing vs. historic no fishing
Mdiff_foref_histp = (Msum_foref - Msum_histp) ./ Msum_histp;
Tdiff_foref_histp = (Tsum_foref - Tsum_histp) ./ Tsum_histp;

%Future fishing vs. historic fishing
Mdiff_foref_histf = (Msum_foref - Msum_histf) ./ Msum_histf;
Tdiff_foref_histf = (Tsum_foref - Tsum_histf) ./ Tsum_histf;

%Future no fishing vs. future no fishing
Mdiff_foref_forep = (Msum_foref - Msum_forep) ./ Msum_forep;
Tdiff_foref_forep = (Tsum_foref - Tsum_forep) ./ Tsum_forep;

%Future no fishing vs. historic no fishing
Mdiff_forep_histp = (Msum_forep - Msum_histp) ./ Msum_histp;
Tdiff_forep_histp = (Tsum_forep - Tsum_histp) ./ Tsum_histp;

% Historic vs. Pre-industrial
hpc(1,:) = Mdiff_histp_pre;
hpc(2,:) = Mdiff_histf_histp;
hpc(3,:) = Mdiff_histf_pre;

% Future vs. Historic 
fhc(1,:) = Mdiff_forep_histp;
fhc(2,:) = Mdiff_foref_forep;
fhc(3,:) = Mdiff_foref_histf;
fhc(4,:) = Mdiff_foref_histp;

%% Bar plots by simulation

%Total
figure(6)
bar(hpc(:,1)*100)
colormap(cmap3)
ylim([-20 0])
set(gca,'XTickLabel',{'Indust CO_2','Fishing','Indust CO_2 + Fishing'})
ylabel('Percent change')
%xlabel('Driver')
title('Global change in mean fish biomass 1951-2000 vs. 1801-1850')
print('-dpng',[ppath 'Global_mean_biom_pre_hist.png'])

%Type
figure(7)
bar(hpc(:,2:4)*100)
colormap(cmap4)
ylim([-50 50])
set(gca,'XTickLabel',{'Indust CO_2','Fishing','Indust CO_2 + Fishing'})
ylabel('Percent change')
legend('Forage','Pelagic','Demersal')
legend('location','northwest')
title('Global change in mean fish biomass 1951-2000 vs. 1801-1850')
print('-dpng',[ppath 'Global_mean_biom_type_pre_hist.png'])

%Size
figure(8)
bar(hpc(:,5:7)*100)
colormap(cmap4)
ylim([-50 50])
set(gca,'XTickLabel',{'Indust CO_2','Fishing','Indust CO_2 + Fishing'})
ylabel('Percent change')
legend('S','M','L')
legend('location','northwest')
title('Global change in mean fish biomass 1951-2000 vs. 1801-1850')
print('-dpng',[ppath 'Global_mean_biom_size_pre_hist.png'])


%Total
figure(9)
bar(fhc(:,1)*100)
colormap(cmap3)
ylim([-20 0])
set(gca,'XTickLabel',{'CC','Fishing','CC + Fishing contemp','CC + Fishing pristine'})
ylabel('Percent change')
%xlabel('Driver')
title('Global change in mean fish biomass 2051-2100 vs. 1951-2000')
print('-dpng',[ppath 'Global_mean_biom_hist_fore.png'])

%Type
figure(10)
bar(fhc(:,2:4)*100)
colormap(cmap4)
ylim([-50 50])
set(gca,'XTickLabel',{'CC','Fishing','CC + Fishing contemp','CC + Fishing pristine'})
ylabel('Percent change')
legend('Forage','Pelagic','Demersal')
legend('location','northwest')
title('Global change in mean fish biomass 2051-2100 vs. 1951-2000')
print('-dpng',[ppath 'Global_mean_biom_type_hist_fore.png'])

%Size
figure(11)
bar(fhc(:,5:7)*100)
colormap(cmap4)
ylim([-50 50])
set(gca,'XTickLabel',{'CC','Fishing','CC + Fishing contemp','CC + Fishing pristine'})
ylabel('Percent change')
legend('S','M','L')
legend('location','northwest')
title('Global change in mean fish biomass 2051-2100 vs. 1951-2000')
print('-dpng',[ppath 'Global_mean_biom_size_hist_fore.png'])

%% Bar plot by size
ML_histf = nansum(all_histf(:,6) + all_histf(:,7));
ML_foref = nansum(all_foref(:,6) + all_foref(:,7));
MLdiff_foref_histf = (ML_foref - ML_histf) ./ ML_histf;

Sdiff(1,:) = -0.036;                  %NPP
Sdiff(2,:) = -0.079;                  %Mesozoop
Sdiff(3,:) = Mdiff_foref_histf(1);    %fish

figure(12)
bar(Sdiff*100)
colormap(cmap3)
%ylim([-40 40])
set(gca,'XTickLabel',{'Net Primary Prod','Mesozooplankton','M&L Fish'})
ylabel('Percent change')
%xlabel('Driver')
title('Global change in productivity')
print('-dpng',[ppath 'Global_prod_size_pre_hist.png'])

%% Time series
load([dpath 'All_biom_mcatch_hist_pristine.mat'],'all_bio');
all_bio_ts = nansum(all_bio);
ts_bio_histp = all_bio_ts(2:16);
clear all_bio all_bio_ts

load([dpath 'All_biom_mcatch_hist_fished.mat'],'all_bio','all_mcatch');
all_bio_ts = nansum(all_bio);
ts_bio_histf = all_bio_ts(2:16);
all_mcatch_ts = nansum(all_mcatch);
ts_mcatch_histf = all_mcatch_ts(2:16);
clear all_bio all_bio_ts all_mcatch all_mcatch_ts

load([dpath 'All_biom_fore_pristine.mat'],'all_bio');
all_bio_ts = nansum(all_bio);
ts_bio_forep = all_bio_ts(2:11);
clear all_bio all_bio_ts

load([dpath 'All_biom_mcatch_fore_fished.mat'],'all_bio','all_mcatch');
all_bio_ts = nansum(all_bio);
ts_bio_foref = all_bio_ts(2:11);
all_mcatch_ts = nansum(all_mcatch);
ts_mcatch_foref = all_mcatch_ts(2:11);
clear all_bio all_bio_ts all_mcatch all_mcatch_ts

yrh=1865:10:2005;
yrf=2005:10:2100;

figure(70)
plot(yrh,ts_bio_histp,'k','LineWidth',2); hold on;
plot(yrh,ts_bio_histf,'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yrf,ts_bio_forep,'k','LineWidth',2); hold on;
plot(yrf,ts_bio_foref,'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
xlim([1860 2100])
xlabel('Year')
ylabel('All fish mean biomass (g/m^2)')
%legend('Pristine','Fished')
%legend('location','southwest')
title('Mean biomass')
print('-dpng',[ppath 'ts_mbio_all.png'])

figure(71)
plot(yrh,ts_mcatch_histf*365,'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yrf,ts_mcatch_foref*365,'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
xlim([1860 2100])
xlabel('Year')
ylabel('All fish mean annual catch (g/m^2)')
%legend('Pristine','Fished')
%legend('location','southwest')
title('Mean catch')
print('-dpng',[ppath 'ts_mcatch_all.png'])

%% N Pac
[npid,ix,grid_np] = intersect(NPid,grid(:,1));

Msum_pre_np = nansum(all_pre(grid_np,:));
Msum_histp_np = nansum(all_histp(grid_np,:));
Msum_histf_np = nansum(all_histf(grid_np,:));
Msum_forep_np = nansum(all_forep(grid_np,:));
Msum_foref_np = nansum(all_foref(grid_np,:));

Tsum_pre_np = nansum(tot_all_pre(grid_np,:));
Tsum_histp_np = nansum(tot_all_histp(grid_np,:));
Tsum_histf_np = nansum(tot_all_histf(grid_np,:));
Tsum_forep_np = nansum(tot_all_forep(grid_np,:));
Tsum_foref_np = nansum(tot_all_foref(grid_np,:));

%% Sims
%Indust vs. Preindust
Mdiff_histp_pre_np = (Msum_histp_np - Msum_pre_np) ./ Msum_pre_np;
Tdiff_histp_pre_np = (Tsum_histp_np - Tsum_pre_np) ./ Tsum_pre_np;
%Hist fishing vs. PreIndust
Mdiff_histf_pre_np = (Msum_histf_np-Msum_pre_np) ./ Msum_pre_np;
Tdiff_histf_pre_np = (Tsum_histf_np-Tsum_pre_np) ./ Tsum_pre_np;
%Hist fishing vs. no fishing
Mdiff_histf_histp_np = (Msum_histf_np - Msum_histp_np) ./ Msum_histp_np;
Tdiff_histf_histp_np = (Tsum_histf_np - Tsum_histp_np) ./ Tsum_histp_np;
%Future fishing vs. historic fishing
Mdiff_foref_histf_np = (Msum_foref_np - Msum_histf_np) ./ Msum_histf_np;
Tdiff_foref_histf_np = (Tsum_foref_np - Tsum_histf_np) ./ Tsum_histf_np;
%Future fishing vs. historic no fishing
Mdiff_foref_histp_np = (Msum_foref_np-Msum_histp_np) ./ Msum_histp_np;
Tdiff_foref_histp_np = (Tsum_foref_np-Tsum_histp_np) ./ Tsum_histp_np;
%Future no fishing vs. historic no fishing
Mdiff_forep_histp_np = (Msum_forep_np-Msum_histp_np) ./ Msum_histp_np;
Tdiff_forep_histp_np = (Tsum_forep_np-Tsum_histp_np) ./ Tsum_histp_np;
%Future fishing vs. no fishing
Mdiff_foref_forep_np = (Msum_foref_np - Msum_forep_np) ./ Msum_forep_np;
Tdiff_foref_forep_np = (Tsum_foref_np - Tsum_forep_np) ./ Tsum_forep_np;

% NPac biomes

%Global biomes mean biomass in g/m2
pre_b(:,1) = nansum(pre_biome(:,1:8),2);
pre_b(:,2) = pre_biome(:,1) + pre_biome(:,4); %all f
pre_b(:,3) = pre_biome(:,2) + pre_biome(:,5) + pre_biome(:,7); %all p
pre_b(:,4) = pre_biome(:,3) + pre_biome(:,6) + pre_biome(:,8); %all d
pre_b(:,5) = nansum(pre_biome(:,1:3),2); %all s
pre_b(:,6) = nansum(pre_biome(:,4:6),2); %all m
pre_b(:,7) = nansum(pre_biome(:,7:8),2); %all l

histp_b(:,1) = nansum(histp_biome(:,1:8),2);
histp_b(:,2) = histp_biome(:,1) + histp_biome(:,4);
histp_b(:,3) = histp_biome(:,2) + histp_biome(:,5) + histp_biome(:,7);
histp_b(:,4) = histp_biome(:,3) + histp_biome(:,6) + histp_biome(:,8);
histp_b(:,5) = nansum(histp_biome(:,1:3),2);
histp_b(:,6) = nansum(histp_biome(:,4:6),2);
histp_b(:,7) = nansum(histp_biome(:,7:8),2);

histf_b(:,1) = nansum(histf_biome(:,1:8),2);
histf_b(:,2) = histf_biome(:,1) + histf_biome(:,4);
histf_b(:,3) = histf_biome(:,2) + histf_biome(:,5) + histf_biome(:,7);
histf_b(:,4) = histf_biome(:,3) + histf_biome(:,6) + histf_biome(:,8);
histf_b(:,5) = nansum(histf_biome(:,1:3),2);
histf_b(:,6) = nansum(histf_biome(:,4:6),2);
histf_b(:,7) = nansum(histf_biome(:,7:8),2);

foref_b(:,1) = nansum(foref_biome(:,1:8),2);
foref_b(:,2) = foref_biome(:,1) + foref_biome(:,4);
foref_b(:,3) = foref_biome(:,2) + foref_biome(:,5) + foref_biome(:,7);
foref_b(:,4) = foref_biome(:,3) + foref_biome(:,6) + foref_biome(:,8);
foref_b(:,5) = nansum(foref_biome(:,1:3),2);
foref_b(:,6) = nansum(foref_biome(:,4:6),2);
foref_b(:,7) = nansum(foref_biome(:,7:8),2);

% NPac biomes grid cells
%Preindust
pb1 = find(biome_pre==1);
NPpb1 = intersect(NPid,pb1);
[nppb1,ia,grid_nppb1] = intersect(NPpb1,grid(:,1));
pb2 = find(biome_pre==2);
NPpb2 = intersect(NPid,pb2);
[nppb2,ia2,grid_nppb2] = intersect(NPpb2,grid(:,1));
pb3 = find(biome_pre==3);
NPpb3 = intersect(NPid,pb3);
[nppb3,ia3,grid_nppb3] = intersect(NPpb3,grid(:,1));
%Hist
hb1 = find(biome_hist==1);
NPhb1 = intersect(NPid,hb1);
[nphb1,ia4,grid_nphb1] = intersect(NPhb1,grid(:,1));
hb2 = find(biome_hist==2);
NPhb2 = intersect(NPid,hb2);
[nphb2,ia5,grid_nphb2] = intersect(NPhb2,grid(:,1));
hb3 = find(biome_hist==3);
NPhb3 = intersect(NPid,hb3);
[nphb3,ia6,grid_nphb3] = intersect(NPhb3,grid(:,1));
%Fore
fb1 = find(biome_fore==1);
NPfb1 = intersect(NPid,fb1);
[npfb1,ia7,grid_npfb1] = intersect(NPfb1,grid(:,1));
fb2 = find(biome_fore==2);
NPfb2 = intersect(NPid,fb2);
[npfb2,ia8,grid_npfb2] = intersect(NPfb2,grid(:,1));
fb3 = find(biome_fore==3);
NPfb3 = intersect(NPid,fb3);
[npfb3,ia9,grid_npfb3] = intersect(NPfb3,grid(:,1));

NPbiome_pre = zeros(360,200);
NPbiome_hist = NPbiome_pre;
NPbiome_fore = NPbiome_pre;
NPbiome_pre(NPpb1) = ones(size(NPpb1));
NPbiome_hist(NPhb1) = ones(size(NPhb1));
NPbiome_fore(NPfb1) = ones(size(NPfb1));
NPbiome_pre(NPpb2) = 2*ones(size(NPpb2));
NPbiome_hist(NPhb2) = 2*ones(size(NPhb2));
NPbiome_fore(NPfb2) = 2*ones(size(NPfb2));
NPbiome_pre(NPpb3) = 3*ones(size(NPpb3));
NPbiome_hist(NPhb3) = 3*ones(size(NPhb3));
NPbiome_fore(NPfb3) = 3*ones(size(NPfb3));


%% calcs
%Only total g to see effect of change in biome size
sum_pre_np1 = nansum(tot_all_pre(grid_nppb1,:));
sum_pre_np2 = nansum(tot_all_pre(grid_nppb2,:));
sum_pre_np3 = nansum(tot_all_pre(grid_nppb3,:));

sum_histp_np1 = nansum(tot_all_histp(grid_nphb1,:));
sum_histp_np2 = nansum(tot_all_histp(grid_nphb2,:));
sum_histp_np3 = nansum(tot_all_histp(grid_nphb3,:));

sum_histf_np1 = nansum(tot_all_histf(grid_nphb1,:));
sum_histf_np2 = nansum(tot_all_histf(grid_nphb2,:));
sum_histf_np3 = nansum(tot_all_histf(grid_nphb3,:));

sum_forep_np1 = nansum(tot_all_forep(grid_npfb1,:));
sum_forep_np2 = nansum(tot_all_forep(grid_npfb2,:));
sum_forep_np3 = nansum(tot_all_forep(grid_npfb3,:));

sum_foref_np1 = nansum(tot_all_foref(grid_npfb1,:));
sum_foref_np2 = nansum(tot_all_foref(grid_npfb2,:));
sum_foref_np3 = nansum(tot_all_foref(grid_npfb3,:));

%Indust vs. Preindust
diff_histp_pre_np(1,:) = (sum_histp_np1 - sum_pre_np1) ./ sum_pre_np1;
diff_histp_pre_np(2,:) = (sum_histp_np2 - sum_pre_np2) ./ sum_pre_np2;
diff_histp_pre_np(3,:) = (sum_histp_np3 - sum_pre_np3) ./ sum_pre_np3;
%Fishing vs. no fishing hist
diff_histf_histp_np(1,:) = (sum_histf_np1 - sum_histp_np1) ./ sum_histp_np1;
diff_histf_histp_np(2,:) = (sum_histf_np2 - sum_histp_np2) ./ sum_histp_np2;
diff_histf_histp_np(3,:) = (sum_histf_np3 - sum_histp_np3) ./ sum_histp_np3;
%Future vs. historic fishing
diff_foref_histf_np(1,:) = (sum_foref_np1 - sum_histf_np1) ./ sum_histf_np1;
diff_foref_histf_np(2,:) = (sum_foref_np2 - sum_histf_np2) ./ sum_histf_np2;
diff_foref_histf_np(3,:) = (sum_foref_np3 - sum_histf_np3) ./ sum_histf_np3;
%Future vs. historic no fishing
diff_forep_histp_np(1,:) = (sum_forep_np1 - sum_histp_np1) ./ sum_histp_np1;
diff_forep_histp_np(2,:) = (sum_forep_np2 - sum_histp_np2) ./ sum_histp_np2;
diff_forep_histp_np(3,:) = (sum_forep_np3 - sum_histp_np3) ./ sum_histp_np3;

diff_foref_forep_np(1,:) = (sum_foref_np1 - sum_forep_np1) ./ sum_forep_np1;
diff_foref_forep_np(2,:) = (sum_foref_np2 - sum_forep_np2) ./ sum_forep_np2;
diff_foref_forep_np(3,:) = (sum_foref_np3 - sum_forep_np3) ./ sum_forep_np3;

%% Mean g/m2 to see effect of just biomass
mean_histf_np1 = nanmean(all_histf(grid_nphb1,:));
mean_histf_np2 = nanmean(all_histf(grid_nphb2,:));
mean_histf_np3 = nanmean(all_histf(grid_nphb3,:));

mean_foref_np1 = nanmean(all_foref(grid_npfb1,:));
mean_foref_np2 = nanmean(all_foref(grid_npfb2,:));
mean_foref_np3 = nanmean(all_foref(grid_npfb3,:));

%Future vs. historic
Bdiff_foref_histf_np(1,:) = (mean_foref_np1 - mean_histf_np1) ./ mean_histf_np1;
Bdiff_foref_histf_np(2,:) = (mean_foref_np2 - mean_histf_np2) ./ mean_histf_np2;
Bdiff_foref_histf_np(3,:) = (mean_foref_np3 - mean_histf_np3) ./ mean_histf_np3;

%% Get change in biome area
area_histf_np(1) = nansum(grid(grid_nphb1,5));
area_histf_np(2) = nansum(grid(grid_nphb2,5));
area_histf_np(3) = nansum(grid(grid_nphb3,5));

area_foref_np(1) = nansum(grid(grid_npfb1,5));
area_foref_np(2) = nansum(grid(grid_npfb2,5));
area_foref_np(3) = nansum(grid(grid_npfb3,5));

%Future vs. historic
Adiff_foref_histf_np = (area_foref_np - area_histf_np) ./ area_histf_np;

%% NPac Historic vs. Pre-industrial
hpc_np(1,:) = Mdiff_histp_pre_np;
hpc_np(2,:) = Mdiff_histf_histp_np;
hpc_np(3,:) = Mdiff_histf_pre_np;

%Total
figure(6)
bar(hpc_np(:,1)*100)
colormap(cmap3)
ylim([-15 0])
set(gca,'XTickLabel',{'Indust CO_2','Fishing','Indust CO_2 + Fishing'})
ylabel('Percent change')
%xlabel('Driver')
title('N Pac change in mean fish biomass 1951-2000 vs. 1801-1850')
print('-dpng',[ppath 'NPac_mean_biom_pre_hist.png'])

%Type
figure(7)
bar(hpc_np(:,2:4)*100)
colormap(cmap4)
ylim([-50 50])
set(gca,'XTickLabel',{'Indust CO_2','Fishing','Indust CO_2 + Fishing'})
ylabel('Percent change')
legend('Forage','Pelagic','Demersal')
legend('location','northwest')
title('N Pac change in mean fish biomass 1951-2000 vs. 1801-1850')
print('-dpng',[ppath 'NPac_mean_biom_type_pre_hist.png'])

%Size
figure(8)
bar(hpc_np(:,5:7)*100)
colormap(cmap4)
ylim([-50 50])
set(gca,'XTickLabel',{'Indust CO_2','Fishing','Indust CO_2 + Fishing'})
ylabel('Percent change')
legend('S','M','L')
legend('location','northwest')
title('N Pac change in mean fish biomass 1951-2000 vs. 1801-1850')
print('-dpng',[ppath 'NPac_mean_biom_size_pre_hist.png'])


%% Future vs. Historic 
fhc_np(1,:) = Mdiff_forep_histp_np;
fhc_np(2,:) = Mdiff_foref_forep_np;
fhc_np(3,:) = Mdiff_foref_histf_np;
fhc_np(4,:) = Mdiff_foref_histp_np;

%Total
figure(9)
bar(fhc_np(:,1)*100)
colormap(cmap3)
ylim([-15 0])
set(gca,'XTickLabel',{'CC','Fishing','CC + Fishing contemp','CC + Fishing pristine'})
ylabel('Percent change')
%xlabel('Driver')
title('N Pac change in mean fish biomass 2006-2100 vs. 1951-2000')
print('-dpng',[ppath 'NPac_mean_biom_hist_fore.png'])

%Type
figure(10)
bar(fhc_np(:,2:4)*100)
colormap(cmap4)
ylim([-50 50])
set(gca,'XTickLabel',{'CC','Fishing','CC + Fishing contemp','CC + Fishing pristine'})
ylabel('Percent change')
legend('Forage','Pelagic','Demersal')
legend('location','northwest')
title('N Pac change in mean fish biomass 2006-2100 vs. 1951-2000')
print('-dpng',[ppath 'NPac_mean_biom_type_hist_fore.png'])

%Size
figure(11)
bar(fhc_np(:,5:7)*100)
colormap(cmap4)
ylim([-50 50])
set(gca,'XTickLabel',{'CC','Fishing','CC + Fishing contemp','CC + Fishing pristine'})
ylabel('Percent change')
legend('S','M','L')
legend('location','northwest')
title('N Pac change in mean fish biomass 2006-2100 vs. 1951-2000')
print('-dpng',[ppath 'NPac_mean_biom_size_hist_fore.png'])



%% N Pac biomes
Tdiff(1,:) = Tdiff_histp_pre;
Tdiff(2,:) = Tdiff_histf_histp;
Tdiff(3,:) = Tdiff_foref_histf;
Tdiff(4,:) = Tdiff_histp_pre_np;
Tdiff(5:7,:) = diff_histp_pre_np;
Tdiff(8,:) = Tdiff_histf_histp_np;
Tdiff(9:11,:) = diff_histf_histp_np;
Tdiff(12,:) = Tdiff_foref_histf_np;
Tdiff(13:15,:) = diff_foref_histf_np;

TdiffA(1,1) = Tdiff(4,1);
TdiffA(1,2) = Tdiff(5,1);
TdiffA(1,3) = Tdiff(6,1);
TdiffA(1,4) = Tdiff(7,1);
TdiffA(2,1) = Tdiff(8,1);
TdiffA(2,2) = Tdiff(9,1);
TdiffA(2,3) = Tdiff(10,1);
TdiffA(2,4) = Tdiff(11,1);
TdiffA(3,1) = Tdiff(12,1);
TdiffA(3,2) = Tdiff(13,1);
TdiffA(3,3) = Tdiff(14,1);
TdiffA(3,4) = Tdiff(15,1);

% hpc_np(1,:) = diff_histp_pre_np;
% hpc_np(2,:) = diff_histf_histp_np;
% hpc_np(3,:) = diff_histf_pre_np;
% fhc_np(1,:) = diff_forep_histp_np;
% fhc_np(2,:) = diff_foref_forep_np;
% fhc_np(3,:) = diff_foref_histf_np;
% fhc_np(4,:) = diff_foref_histp_np;

figure(36)
bar(TdiffA*100)
colormap(cmap3)
set(gca,'XTickLabel',{'Industrial CO_2','Fishing','Climate Change + Fishing'})
ylabel('Percent change')
xlim([2.5 3.5])
ylim([-30 30])
legend('All','LC','ECCS','ECSS')
legend('location','southwest')
title('N Pac change in mean fish total abundance (g)')
print('-dpng',[ppath 'NPacBiome_tot_mean_biom_foref_histf.png'])

biome_diff(1,1:3)=Adiff_foref_histf_np;
biome_diff(2,1:3)=Bdiff_foref_histf_np(:,1)';

figure(37)
bar(biome_diff*100)
colormap(cmap2)
ylim([-30 30])
set(gca,'XTickLabel',{'Biome area (m)','Fish biomass (g/m^2)'})
ylabel('Percent change')
legend('LC','ECCS','ECSS')
legend('location','southwest')
title('N Pac biome changes')
print('-dpng',[ppath 'NPacBiome_changes_foref_histf.png'])

figure(38)
bar(biome_diff*100)
colormap(cmap2)
xlim([0.5 1.5])
ylim([-30 30])
set(gca,'XTickLabel',{'Biome area (m)','Fish biomass (g/m^2)'})
ylabel('Percent change')
legend('LC','ECCS','ECSS')
legend('location','southwest')
title('N Pac change in biome area')
print('-dpng',[ppath 'NPacBiome_changes_area_foref_histf.png'])


figure(39)
bar(biome_diff*100)
colormap(cmap2)
xlim([1.5 2.5])
ylim([-30 30])
set(gca,'XTickLabel',{'Biome area (m)','Fish biomass (g/m^2)'})
ylabel('Percent change')
legend('LC','ECCS','ECSS')
legend('location','southeast')
title('N Pac biome change in fish mean abundance')
print('-dpng',[ppath 'NPacBiome_changes_biomass_foref_histf.png'])


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

%MOM grid
All_histf_mg = NaN*ones(size(geolon_t));
All_foref_mg = NaN*ones(size(geolon_t));

All_histf_mg(grid(:,1)) = all_histf(:,1);
All_foref_mg(grid(:,1)) = all_foref(:,1);

Diff_foref_histf_gg = (All_foref_gg - All_histf_gg) ./All_histf_gg;
Diff_foref_histf_mg = (All_foref_mg - All_histf_mg) ./All_histf_mg;

%%
figure(21)
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
print('-dpng',[ppath 'Diff_global_all_foref_histf_v1.png'])

%%
figure(22)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,Diff_foref_histf_mg*100)
colormap(cmap_color_rb)              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 50]);
hcb = colorbar('h');
ylim(hcb,[-50 50])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference from 1951-2000 mean biomass of all fishes')
axesmui
stamp(cfile)
print('-dpng',[ppath 'Diff_global_all_foref_histf_caxis.png'])

%% all L
figure(23)
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

%% Biomes
figure(50)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,NPbiome_pre)
colormap(cmap1)              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gcf,'renderer','painters')
title('N Pac Biomes Pre-Industrial 1801-1850')
print('-dpng',[cpath 'NPac_biomes_preindust.png'])

figure(51)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,NPbiome_hist)
colormap(cmap1)          
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gcf,'renderer','painters')
title('N Pac Biomes Historic 1951-2000')
print('-dpng',[cpath 'NPac_biomes_historic.png'])

figure(52)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,NPbiome_fore)
colormap(cmap1)          
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gcf,'renderer','painters')
title('N Pac Biomes Forecast 2051-2100')
print('-dpng',[cpath 'NPac_biomes_forecast.png'])

%% Size & type Global and NPac together
% Future vs. Historic 
fhc_both(1,:) = Mdiff_foref_histf;
fhc_both(2,:) = Mdiff_foref_histf_np;

%Type
figure(10)
bar(fhc_both(:,2:4)*100)
colormap(cmap4)
ylim([-30 10])
set(gca,'XTickLabel',{'Global','N Pac'})
ylabel('Percent change')
legend('Forage','Pelagic','Demersal')
legend('location','northwest')
title('Change in mean fish biomass 2006-2100 vs. 1951-2000')
print('-dpng',[ppath 'both_mean_biom_type_hist_fore.png'])

%Size
figure(11)
bar(fhc_both(:,5:7)*100)
colormap(cmap4)
ylim([-30 10])
set(gca,'XTickLabel',{'Global','N Pac'})
ylabel('Percent change')
legend('S','M','L')
legend('location','northwest')
title('Change in mean fish biomass 2006-2100 vs. 1951-2000')
print('-dpng',[ppath 'both_mean_biom_size_hist_fore.png'])

figure(12)
subplot(1,2,1)
bar(fhc_both(:,2:4)*100)
colormap(cmap4)
ylim([-30 10])
set(gca,'XTickLabel',{'Global','N Pac'})
ylabel('Percent change')
legend('Forage','Pelagic','Demersal')
legend('location','northwest')
text(2.75,12,'Change in mean fish biomass 2006-2100 vs. 1951-2000',...
    'HorizontalAlignment','center','FontWeight','bold')
subplot(1,2,2)
bar(fhc_both(:,5:7)*100)
colormap(cmap4)
ylim([-30 10])
set(gca,'XTickLabel',{'Global','N Pac'})
ylabel('Percent change')
legend('S','M','L')
legend('location','northwest')
print('-dpng',[ppath 'both_mean_biom_size_type_hist_fore.png'])

%% By biome type & size

%Type
figure(13)
bar(diff_foref_histf_np(:,2:4)'*100)
colormap(cmap2)
ylim([-50 40])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
legend('LC','ECCS','ECSS')
legend('location','northwest')
title('Change in mean fish biomass 2006-2100 vs. 1951-2000')
print('-dpng',[ppath 'NPacBiome_biom_type_hist_fore.png'])

%Size
figure(14)
bar(diff_foref_histf_np(:,5:7)'*100)
colormap(cmap2)
ylim([-50 40])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
legend('LC','ECCS','ECSS')
legend('location','northwest')
title('Change in mean fish biomass 2006-2100 vs. 1951-2000')
print('-dpng',[ppath 'NPacBiome_biom_size_hist_fore.png'])
%%
figure(15)
subplot(1,2,1)
bar(diff_foref_histf_np(:,2:4)'*100)
colormap(cmap2)
ylim([-50 40])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
legend('LC','ECCS','ECSS')
legend('location','northwest')
text(4.5,45,'Change in mean fish biomass 2006-2100 vs. 1951-2000',...
    'HorizontalAlignment','center','FontWeight','bold')
subplot(1,2,2)
bar(diff_foref_histf_np(:,5:7)'*100)
colormap(cmap2)
ylim([-50 40])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
print('-dpng',[ppath 'NPacBiome_biom_size_type_hist_fore.png'])


