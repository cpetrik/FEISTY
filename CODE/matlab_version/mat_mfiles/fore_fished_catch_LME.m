% Changes to catch by LMEs 

clear all
close all

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
dp = '/Volumes/GFDL/NC/Jul_og_sizes/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Julia_OG_sizes/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([gpath 'hindcast_gridspec.mat'],'dat','geolat_t','geolon_t');
load([gpath 'lme_mask_esm2m.mat']);
load([gpath 'NPid_esm2m.mat']);
grid = csvread([gpath 'grid_csv.csv']);

%Colormap
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

%% Hist fished
load([dpath 'LME_hist_fished_' cfile '.mat'],'lme_mcatch00');
load([dpath 'Means_hist_fished_' cfile '.mat'],'mf_mcatch5000','mp_mcatch5000',...
    'md_mcatch5000','lp_mcatch5000','ld_mcatch5000');

histf_lme = lme_mcatch00;

%Mean global biomass in g/m2
all_histf(:,1) = mf_mcatch5000+...
    mp_mcatch5000+md_mcatch5000+lp_mcatch5000+ld_mcatch5000;
all_histf(:,2) = mf_mcatch5000; %all f
all_histf(:,3) = mp_mcatch5000+lp_mcatch5000; %all p
all_histf(:,4) = md_mcatch5000+ld_mcatch5000; %all d

%Mean total global catch in g
tot_all_histf = all_histf .* repmat(grid(:,5),1,4);

clear lme_mcatch00 sf_mcatch5000 sp_mcatch5000 sd_mcatch5000
clear mf_mcatch5000 mp_mcatch5000 md_mcatch5000 b_mcatch5000 lp_mcatch5000 ld_mcatch5000

%% Fore fished
load([dpath 'LME_fore_fished_' cfile '.mat'],'lme_mcatch00');
load([dpath 'Means_fore_fished_' cfile '.mat'],'mf_mcatch5000','mp_mcatch5000',...
    'md_mcatch5000','lp_mcatch5000','ld_mcatch5000');

foref_lme = lme_mcatch00;

%Mean global biomass in g/m2
all_foref(:,1) = mf_mcatch5000+...
    mp_mcatch5000+md_mcatch5000+lp_mcatch5000+ld_mcatch5000;
all_foref(:,2) = mf_mcatch5000; %all f
all_foref(:,3) = mp_mcatch5000+lp_mcatch5000; %all p
all_foref(:,4) = md_mcatch5000+ld_mcatch5000; %all d

%Mean total global catch in g
tot_all_foref = all_foref .* repmat(grid(:,5),1,4);

clear lme_mcatch00 sf_mcatch5000 sp_mcatch5000 sd_mcatch5000
clear mf_mcatch5000 mp_mcatch5000 md_mcatch5000 b_mcatch5000 lp_mcatch5000 ld_mcatch5000

%% Global
%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;
x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);

Call_hist=griddata(grid(:,2),grid(:,3),all_histf(:,1),X,Y);
CF_hist=griddata(grid(:,2),grid(:,3),all_histf(:,2),X,Y);
CP_hist=griddata(grid(:,2),grid(:,3),all_histf(:,3),X,Y);
CD_hist=griddata(grid(:,2),grid(:,3),all_histf(:,4),X,Y);

Call_fore=griddata(grid(:,2),grid(:,3),all_foref(:,1),X,Y);
CF_fore=griddata(grid(:,2),grid(:,3),all_foref(:,2),X,Y);
CP_fore=griddata(grid(:,2),grid(:,3),all_foref(:,3),X,Y);
CD_fore=griddata(grid(:,2),grid(:,3),all_foref(:,4),X,Y);

%Future vs. historic
Cdiff_foref_histf = (Call_fore-Call_hist) ./ Call_hist;
FCdiff_foref_histf = (CF_fore-CF_hist) ./ CF_hist;
PCdiff_foref_histf = (CP_fore-CP_hist) ./ CP_hist;
DCdiff_foref_histf = (CD_fore-CD_hist) ./ CD_hist;

%% N Pac
[npid,ix,grid_np] = intersect(NPid,grid(:,1));

Msum_histf_np = nansum(all_histf(grid_np,:));
Msum_foref_np = nansum(all_foref(grid_np,:));

Tsum_histf_np = nansum(tot_all_histf(grid_np,:));
Tsum_foref_np = nansum(tot_all_foref(grid_np,:));

%Future vs. historic
Mdiff_foref_histf_np = (Msum_foref_np - Msum_histf_np) ./ Msum_histf_np;
Tdiff_foref_histf_np = (Tsum_foref_np - Tsum_histf_np) ./ Tsum_histf_np;


%% LMEs
tlme = lme_mask_esm2m';

hlme_mf = NaN*ones(360,200);
hlme_mp = hlme_mf;
hlme_md = hlme_mf;
hlme_lp = hlme_mf;
hlme_ld = hlme_mf;
flme_mf = hlme_mf;
flme_mp = hlme_mf;
flme_md = hlme_mf;
flme_lp = hlme_mf;
flme_ld = hlme_mf;

for L=1:66
    lid = find(tlme==L);
    
    hlme_mf(lid) = histf_lme(L,1);
    hlme_mp(lid) = histf_lme(L,2);
    hlme_md(lid) = histf_lme(L,3);
    hlme_lp(lid) = histf_lme(L,4);
    hlme_ld(lid) = histf_lme(L,5);
    
    flme_mf(lid) = foref_lme(L,1);
    flme_mp(lid) = foref_lme(L,2);
    flme_md(lid) = foref_lme(L,3);
    flme_lp(lid) = foref_lme(L,4);
    flme_ld(lid) = foref_lme(L,5);
end

hlme_All = hlme_mf+hlme_mp+hlme_md+hlme_lp+hlme_ld;
hlme_AllF = hlme_mf;
hlme_AllP = hlme_mp+hlme_lp;
hlme_AllD = hlme_md+hlme_ld;

flme_All = flme_mf+flme_mp+flme_md+flme_lp+flme_ld;
flme_AllF = flme_mf;
flme_AllP = flme_mp+flme_lp;
flme_AllD = flme_md+flme_ld;

Diff_foref_histf = (flme_All-hlme_All) ./ hlme_All;
FDiff_foref_histf = (flme_AllF-hlme_AllF) ./ hlme_AllF;
PDiff_foref_histf = (flme_AllP-hlme_AllP) ./ hlme_AllP;
DDiff_foref_histf = (flme_AllD-hlme_AllD) ./ hlme_AllD;

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

%% Global LME
% ALL
cjrev=colormap('jet');
cjrev=flipud(cjrev);

figure(1)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,Diff_foref_histf*100)
colormap(cjrev)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 30]);
hcb = colorbar('h');
ylim(hcb,[-50 30])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference in catch of all fishes from 1951-2000')
stamp(cfile)
print('-dpng',[ppath 'Diff_global_lme_catch_all_foref_histf.png'])

% all F
figure(2)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,FDiff_foref_histf*100)
colormap(cjrev)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 30]);
hcb = colorbar('h');
ylim(hcb,[-50 30])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference in forage fish catch from 1951-2000')
stamp(cfile)
print('-dpng',[ppath 'Diff_global_lme_catch_F_foref_histf.png'])

% all P
figure(3)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,PDiff_foref_histf*100)
colormap(cjrev)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 30]);
hcb = colorbar('h');
ylim(hcb,[-50 30])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference in large pelagics catch from 1951-2000')
stamp(cfile)
print('-dpng',[ppath 'Diff_global_lme_catch_P_foref_histf.png'])

% All D
figure(4)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,DDiff_foref_histf*100)
colormap(cjrev)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 30]);
hcb = colorbar('h');
ylim(hcb,[-50 30])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference in demersal catch from 1951-2000')
stamp(cfile)
print('-dpng',[ppath 'Diff_global_lme_catch_D_foref_histf.png'])


%% NPac only LME
% ALL
figure(5)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,Diff_foref_histf*100)
colormap(cjrev)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 30]);
hcb = colorbar('h');
ylim(hcb,[-50 30])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference in catch of all fishes from 1951-2000')
stamp(cfile)
print('-dpng',[ppath 'Diff_NPac_lme_catch_all_foref_histf.png'])

% all F
figure(6)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,FDiff_foref_histf*100)
colormap(cjrev)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 30]);
hcb = colorbar('h');
ylim(hcb,[-50 30])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference in forage fish catch from 1951-2000')
stamp(cfile)
print('-dpng',[ppath 'Diff_NPac_lme_catch_F_foref_histf.png'])

% all P
figure(7)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,PDiff_foref_histf*100)
colormap(cjrev)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 30]);
hcb = colorbar('h');
ylim(hcb,[-50 30])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference in large pelagic fish catch from 1951-2000')
stamp(cfile)
print('-dpng',[ppath 'Diff_NPac_lme_catch_P_foref_histf.png'])

% All D
figure(8)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,DDiff_foref_histf*100)
colormap(cjrev)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-50 30]);
hcb = colorbar('h');
ylim(hcb,[-50 30])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference in demersal catch from 1951-2000')
stamp(cfile)
print('-dpng',[ppath 'Diff_NPac_lme_catch_D_foref_histf.png'])

%% Global
figure(9)
m_proj('miller','lat',82);
m_pcolor(X,Y,Cdiff_foref_histf*100); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2051-2100 % difference in catch of all fishes from 1951-2000')
colormap(cjrev)
colorbar('h')
caxis([-100 50])
stamp(cfile)
print('-dpng',[ppath 'Diff_global_catch_all_foref_histf.png'])

figure(10)
m_proj('miller','lat',82);
m_pcolor(X,Y,FCdiff_foref_histf*100); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2051-2100 % difference in forage fish catch from 1951-2000')
colormap(cjrev)
colorbar('h')
caxis([-100 50])
stamp(cfile)
print('-dpng',[ppath 'Diff_global_catch_F_foref_histf.png'])

figure(11)
m_proj('miller','lat',82);
m_pcolor(X,Y,PCdiff_foref_histf*100); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2051-2100 % difference in large pelagics catch from 1951-2000')
colormap(cjrev)
colorbar('h')
caxis([-100 50])
stamp(cfile)
print('-dpng',[ppath 'Diff_global_catch_P_foref_histf.png'])

figure(12)
m_proj('miller','lat',82);
m_pcolor(X,Y,DCdiff_foref_histf*100); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2051-2100 % difference in demersal catch from 1951-2000')
colormap(cjrev)
colorbar('h')
caxis([-100 50])
stamp(cfile)
print('-dpng',[ppath 'Diff_global_catch_D_foref_histf.png'])


