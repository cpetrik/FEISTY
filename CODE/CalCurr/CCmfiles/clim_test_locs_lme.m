% Visualize output of POEM Climatology at single locations
% 150 years, monthly means saved

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
gpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
cpath = ['/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/'];
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([gpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'J&C15_all.mat'],'LMEname')

%LMEs of interest
%[GB','WSS','CSS','ESS'],'GS','NS','NwS','BS','EBS','K2','S1','HOT','BATS','PUp'
locs=[8;19;22;21;20;1;51;49;10;6;13];
sites = LMEname(locs);
name = {'SS','GS','NS','NwS','BS','EBS','OyaCur','KurCur','Haw','SEUS','Humb'};

% Colors
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
cmap3=cmap_ppt([3,1,5],:);
cm={[1 0.5 0],...   %orange
    [0.5 0.5 0],... %tan/army
    [0 0.7 0],...   %g
    [0 1 1],...     %c
    [0 0 0.75],...  %b
    [0.5 0 1],...   %purple
    [1 0 1],...     %m
    [1 0 0],...     %r
    [0.5 0 0],...   %maroon
    [0.75 0.75 0.75],... %lt grey
    [0.5 0.5 0.5],...    %med grey
    [49/255 79/255 79/255],... %dk grey
    [0 0 0],...      %black
    [1 1 0],...      %yellow
    [127/255 255/255 0],... %lime green
    [0 0.5 0],...    %dk green
    [0/255 206/255 209/255],... %turq
    [0 0.5 0.75],...   %med blue
    [188/255 143/255 143/255],... %rosy brown
    [255/255 192/255 203/255],... %pink
    [255/255 160/255 122/255]}; %peach

% plot info
[ni,nj]=size(lon);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

geolat_t=lat;
geolon_t=lon;

%% POEM means
cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';
ppath = [pp cfile '/'];
dpath = [dp cfile '/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat'])
load([dpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);
load([dpath 'Means_con_types_Climatol_' harv '_' cfile '.mat']);


%% Zoop, det, bent
load([cpath 'cobalt_zoop_biom_means.mat'],'mz_mean_clim','lz_mean_clim','mzloss_mean_clim','lzloss_mean_clim')
load([cpath 'cobalt_det_biom_means.mat'],'det_mean_clim')

%ESM2.6 in mg C m-2 or mg C m-2 d-1
%from mg C m-2 to g(WW) m-2
% 1e-3 g C in 1 mg C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
zmean_grid = (mz_mean_clim + lz_mean_clim) * 1e-3 * 9.0;
zloss_grid = (mzloss_mean_clim+lzloss_mean_clim) * 1e-3 * 9.0;
det_grid = det_mean_clim * 1e-3 * 9.0;

%% Grid
Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zb=NaN*ones(ni,nj);

%mean
Zsf(ID)=sf_mean;
Zsp(ID)=sp_mean;
Zsd(ID)=sd_mean;
Zmf(ID)=mf_mean;
Zmp(ID)=mp_mean;
Zmd(ID)=md_mean;
Zlp(ID)=lp_mean;
Zld(ID)=ld_mean;
Zb(ID)=b_mean;

Csf=NaN*ones(ni,nj);
Csp=NaN*ones(ni,nj);
Csd=NaN*ones(ni,nj);
CmfZm=NaN*ones(ni,nj);
CmfZl=NaN*ones(ni,nj);
CmfF=NaN*ones(ni,nj);
CmfP=NaN*ones(ni,nj);
CmpZm=NaN*ones(ni,nj);
CmpZl=NaN*ones(ni,nj);
CmpF=NaN*ones(ni,nj);
CmpP=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
ClpF=NaN*ones(ni,nj);
ClpP=NaN*ones(ni,nj);
CldF=NaN*ones(ni,nj);
CldP=NaN*ones(ni,nj);
CldD=NaN*ones(ni,nj);
CldB=NaN*ones(ni,nj);

%
Csf(ID)=sf_mconZm;
Csp(ID)=sp_mconZm;
Csd(ID)=sd_mconZm;

CmfZm(ID)=mf_mconZm;
CmfZl(ID)=mf_mconZl;
CmfF(ID)=mf_mconF;
CmfP(ID)=mf_mconP;

CmpZm(ID)=mp_mconZm;
CmpZl(ID)=mp_mconZl;
CmpF(ID)=mp_mconF;
CmpP(ID)=mp_mconP;

Cmd(ID)=md_mconB;

ClpF(ID)=lp_mconF;
ClpP(ID)=lp_mconP;

CldF(ID)=ld_mconF;
CldP(ID)=ld_mconP;
CldD(ID)=ld_mconD;
CldB(ID)=ld_mconB;

All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllD = Zsd+Zmd+Zld;

CAll = Csf+CmfZm+CmfZl+CmfF+CmfP+...
    Csp+CmpZm+CmpZl+CmpF+CmpP+ClpF+ClpP+...
    Csd+Cmd+CldF+CldP+CldD+CldB;
CAllF = Csf+CmfZm+CmfZl+CmfF+CmfP;
CAllP = Csp+CmpZm+CmpZl+CmpF+CmpP+ClpF+ClpP;
CAllD = Csd+Cmd+CldF+CldP+CldD+CldB;
CAllS = Csf+Csp+Csd;
CAllM = CmfZm+CmfZl+CmfF+CmfP+CmpZm+CmpZl+CmpF+CmpP+Cmd;
CAllL = ClpF+ClpP+CldF+CldP+CldD+CldB;


%% Calc LMEs
tlme = lme_mask_onedeg;

lme_biom = NaN*ones(66,4);
lme_prey = NaN*ones(66,4);
lme_cF = zeros(66,6);
lme_cP = zeros(66,6);
lme_cD = zeros(66,6);
for L=1:66
    lid = find(tlme==L);
    %prey
    lme_prey(L,1) = nanmean(zmean_grid(lid));
    lme_prey(L,2) = nanmean(zloss_grid(lid));
    lme_prey(L,3) = nanmean(det_grid(lid));
    lme_prey(L,4) = nanmean(Zb(lid));
    %biomass
    lme_biom(L,1) = nanmean(AllF(lid));
    lme_biom(L,2) = nanmean(AllP(lid));
    lme_biom(L,3) = nanmean(AllD(lid));
    lme_biom(L,4) = nanmean(All(lid));
    %consumption
    lme_cF(L,1) = nanmean(Csf(lid)) + nanmean(CmfZm(lid));
    lme_cF(L,2) = nanmean(CmfZl(lid));
    lme_cF(L,3) = nanmean(CmfF(lid));
    lme_cF(L,4) = nanmean(CmfP(lid));
    
    lme_cP(L,1) = nanmean(Csp(lid)) + nanmean(CmpZm(lid));
    lme_cP(L,2) = nanmean(CmpZl(lid));
    lme_cP(L,3) = nanmean(CmpF(lid));
    lme_cP(L,4) = nanmean(CmpP(lid));
    
    lme_cD(L,1) = nanmean(Csd(lid));
    lme_cD(L,3) = nanmean(CldF(lid));
    lme_cD(L,4) = nanmean(CldP(lid));
    lme_cD(L,5) = nanmean(CldD(lid));
    lme_cD(L,6) = nanmean(Cmd(lid)) + nanmean(CldB(lid));
    
end


%% Table
Tab=table(locs,sites,lme_prey(locs,1),lme_prey(locs,2),lme_prey(locs,3),...
    lme_prey(locs,4),lme_biom(locs,1),lme_biom(locs,2),lme_biom(locs,3),...
    'VariableNames',{'LME','Name','Z','Zloss','Det','Bent','F','P','D'});
writetable(Tab,[dpath 'LME_locs_biomasses_clim_fished_',harv,'_' cfile '.csv'],'Delimiter',',');
save([dpath 'LME_locs_biomasses_clim_fished_',harv,'_' cfile '.mat'],'Tab');

%writetable(Tab,['LME_locs_biomasses_clim_fished_',harv,'_' cfile '.csv'],'Delimiter',',');


%% Plots
x1=linspace(1.2,2.8,9);
x2=linspace(3.2,4.8,9);
x3=linspace(1.2,4.8,9);
x4=5*ones(9,1);
y1=ones(9,1);
y2=2*ones(9,1);
y3=linspace(1.8,1.2,9);

for n=1:length(sites)
    s = locs(n);
    loc = name{n};
    
    %% Bubble/arrow plot ----------------------------------------------------
    
    figure(1)
    clf
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',10*lme_prey(s,1)); hold on;
    plot(1,1,'.','color','k','MarkerSize',10*lme_prey(s,3)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',10*lme_biom(s,1)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',10*lme_biom(s,2)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',10*lme_biom(s,3)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',50*(lme_cF(s,1)+lme_cF(s,2))); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',50*lme_cP(s,3)); hold on;
    plot(x3,y1,'color','k','Linewidth',lme_cD(s,6)); hold on;
    if (lme_cD(s,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',10*lme_cD(s,3)); hold on;
    end
    if (lme_cD(s,4) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',10*lme_cD(s,4)); hold on;
    end
    text(1,2.4,sprintf('%2.1f',lme_prey(s,1))); hold on;
    text(1,0.6,sprintf('%2.1f',lme_prey(s,3))); hold on;
    text(3,2.4,sprintf('%2.1f',lme_biom(s,1))); hold on;
    text(5,2.4,sprintf('%2.1f',lme_biom(s,2))); hold on;
    text(5,0.6,sprintf('%2.1f',lme_biom(s,3))); hold on;
    text(2,2.2,sprintf('%2.1e',(lme_cF(s,1)+lme_cF(s,2)))); hold on;
    text(4,2.2,sprintf('%2.1e',lme_cP(s,3))); hold on;
    text(3,0.8,sprintf('%2.1e',lme_cD(s,6))); hold on;
    text(4,1.5,sprintf('%2.1e',lme_cD(s,3))); hold on;
    text(5.25,1.5,sprintf('%2.1e',lme_cD(s,4))); hold on;
    xlim([0 6])
    ylim([0 3])
    title(LMEname{s})
    set(gca,'XTickLabel','','YTickLabel','')
    print('-dpng',[ppath 'Clim_fished_',harv,'_fluxes_Zbio_det_' loc '_LME.png'])
    
    figure(2)
    clf
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',10*lme_prey(s,2)); hold on;
    plot(1,1,'.','color','k','MarkerSize',10*lme_prey(s,3)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',10*lme_biom(s,1)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',10*lme_biom(s,2)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',10*lme_biom(s,3)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',50*(lme_cF(s,1)+lme_cF(s,2))); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',50*lme_cP(s,3)); hold on;
    plot(x3,y1,'color','k','Linewidth',lme_cD(s,6)); hold on;
    if (lme_cD(s,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',10*lme_cD(s,3)); hold on;
    end
    if (lme_cD(s,4) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',10*lme_cD(s,4)); hold on;
    end
    text(1,2.4,sprintf('%2.1f',lme_prey(s,2))); hold on;
    text(1,0.6,sprintf('%2.1f',lme_prey(s,3))); hold on;
    text(3,2.4,sprintf('%2.1f',lme_biom(s,1))); hold on;
    text(5,2.4,sprintf('%2.1f',lme_biom(s,2))); hold on;
    text(5,0.6,sprintf('%2.1f',lme_biom(s,3))); hold on;
    text(2,2.2,sprintf('%2.1e',(lme_cF(s,1)+lme_cF(s,2)))); hold on;
    text(4,2.2,sprintf('%2.1e',lme_cP(s,3))); hold on;
    text(3,0.8,sprintf('%2.1e',lme_cD(s,6))); hold on;
    text(4,1.5,sprintf('%2.1e',lme_cD(s,3))); hold on;
    text(5.25,1.5,sprintf('%2.1e',lme_cD(s,4))); hold on;
    xlim([0 6])
    ylim([0 3])
    title(LMEname{s})
    set(gca,'XTickLabel','','YTickLabel','')
    print('-dpng',[ppath 'Clim_fished_',harv,'_fluxes_Zloss_det_' loc '_LME.png'])
    
    figure(3)
    clf
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',10*lme_prey(s,1)); hold on;
    plot(1,1,'.','color','k','MarkerSize',10*lme_prey(s,4)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',10*lme_biom(s,1)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',10*lme_biom(s,2)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',10*lme_biom(s,3)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',50*(lme_cF(s,1)+lme_cF(s,2))); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',50*lme_cP(s,3)); hold on;
    plot(x3,y1,'color','k','Linewidth',lme_cD(s,6)); hold on;
    if (lme_cD(s,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',10*lme_cD(s,3)); hold on;
    end
    if (lme_cD(s,4) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',10*lme_cD(s,4)); hold on;
    end
    text(1,2.4,sprintf('%2.1f',lme_prey(s,1))); hold on;
    text(1,0.6,sprintf('%2.1f',lme_prey(s,4))); hold on;
    text(3,2.4,sprintf('%2.1f',lme_biom(s,1))); hold on;
    text(5,2.4,sprintf('%2.1f',lme_biom(s,2))); hold on;
    text(5,0.6,sprintf('%2.1f',lme_biom(s,3))); hold on;
    text(2,2.2,sprintf('%2.1e',(lme_cF(s,1)+lme_cF(s,2)))); hold on;
    text(4,2.2,sprintf('%2.1e',lme_cP(s,3))); hold on;
    text(3,0.8,sprintf('%2.1e',lme_cD(s,6))); hold on;
    text(4,1.5,sprintf('%2.1e',lme_cD(s,3))); hold on;
    text(5.25,1.5,sprintf('%2.1e',lme_cD(s,4))); hold on;
    xlim([0 6])
    ylim([0 3])
    title(LMEname{s})
    set(gca,'XTickLabel','','YTickLabel','')
    print('-dpng',[ppath 'Clim_fished_',harv,'_fluxes_Zbio_bent_' loc '_LME.png'])
    
    figure(4)
    clf
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',10*lme_prey(s,2)); hold on;
    plot(1,1,'.','color','k','MarkerSize',10*lme_prey(s,4)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',10*lme_biom(s,1)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',10*lme_biom(s,2)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',10*lme_biom(s,3)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',50*(lme_cF(s,1)+lme_cF(s,2))); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',50*lme_cP(s,3)); hold on;
    plot(x3,y1,'color','k','Linewidth',lme_cD(s,6)); hold on;
    if (lme_cD(s,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',10*lme_cD(s,3)); hold on;
    end
    if (lme_cD(s,4) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',10*lme_cD(s,4)); hold on;
    end
    text(1,2.4,sprintf('%2.1f',lme_prey(s,2))); hold on;
    text(1,0.6,sprintf('%2.1f',lme_prey(s,4))); hold on;
    text(3,2.4,sprintf('%2.1f',lme_biom(s,1))); hold on;
    text(5,2.4,sprintf('%2.1f',lme_biom(s,2))); hold on;
    text(5,0.6,sprintf('%2.1f',lme_biom(s,3))); hold on;
    text(2,2.2,sprintf('%2.1e',(lme_cF(s,1)+lme_cF(s,2)))); hold on;
    text(4,2.2,sprintf('%2.1e',lme_cP(s,3))); hold on;
    text(3,0.8,sprintf('%2.1e',lme_cD(s,6))); hold on;
    text(4,1.5,sprintf('%2.1e',lme_cD(s,3))); hold on;
    text(5.25,1.5,sprintf('%2.1e',lme_cD(s,4))); hold on;
    xlim([0 6])
    ylim([0 3])
    title(LMEname{s})
    set(gca,'XTickLabel','','YTickLabel','')
    print('-dpng',[ppath 'Clim_fished_',harv,'_fluxes_Zloss_bent_' loc '_LME.png'])
    
end

%% Subplots of each type
%locs=[8;19;22;21;20;1;51;49;10;6;13];
shelf = [8;19;22;21;20;1]; %1-6
upw = 13; %11
olig = [49;10;6];%8-10
fave = [22;10;13];

% Shelf Sea
for n=1:6
    s = locs(n);
    loc = name{n};
    
    f5=figure(5);
    subplot(3,2,n)
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',10*lme_prey(s,1)); hold on;
    plot(1,1,'.','color','k','MarkerSize',10*lme_prey(s,4)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',10*lme_biom(s,1)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',10*lme_biom(s,2)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',10*lme_biom(s,3)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',50*(lme_cF(s,1)+lme_cF(s,2))); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',50*lme_cP(s,3)); hold on;
    plot(x3,y1,'color','k','Linewidth',lme_cD(s,6)); hold on;
    if (lme_cD(s,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',10*lme_cD(s,3)); hold on;
    end
    if (lme_cD(s,4) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',10*lme_cD(s,4)); hold on;
    end
    text(1,2.4,sprintf('%2.1f',lme_prey(s,1))); hold on;
    text(1,0.6,sprintf('%2.1f',lme_prey(s,4))); hold on;
    text(3,2.4,sprintf('%2.1f',lme_biom(s,1))); hold on;
    text(5,2.4,sprintf('%2.1f',lme_biom(s,2))); hold on;
    text(5,0.6,sprintf('%2.1f',lme_biom(s,3))); hold on;
    text(2,2.2,sprintf('%2.1e',(lme_cF(s,1)+lme_cF(s,2)))); hold on;
    text(4,2.2,sprintf('%2.1e',lme_cP(s,3))); hold on;
    text(3,0.8,sprintf('%2.1e',lme_cD(s,6))); hold on;
    text(4,1.5,sprintf('%2.1e',lme_cD(s,3))); hold on;
    text(5.25,1.5,sprintf('%2.1e',lme_cD(s,4))); hold on;
    xlim([0 6])
    ylim([0 3])
    title(LMEname{s})
    set(gca,'XTickLabel','','YTickLabel','')
    
end
print(f5,'-dpng',[ppath 'Clim_fished_',harv,'_fluxes_Zbio_bent_shelf_LMEs.png'])


%% Oligotrophic
count=0;
for n=8:10
    s = locs(n);
    loc = name{n};
    count = count +1;
    
    f6=figure(6);
    subplot(2,2,count)
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',10*lme_prey(s,1)); hold on;
    plot(1,1,'.','color','k','MarkerSize',10*lme_prey(s,4)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',10*lme_biom(s,1)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',10*lme_biom(s,2)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',10*lme_biom(s,3)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',50*(lme_cF(s,1)+lme_cF(s,2))); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',50*lme_cP(s,3)); hold on;
    plot(x3,y1,'color','k','Linewidth',lme_cD(s,6)); hold on;
    if (lme_cD(s,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',10*lme_cD(s,3)); hold on;
    end
    if (lme_cD(s,4) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',10*lme_cD(s,4)); hold on;
    end
    text(1,2.4,sprintf('%2.1f',lme_prey(s,1))); hold on;
    text(1,0.6,sprintf('%2.1f',lme_prey(s,4))); hold on;
    text(3,2.4,sprintf('%2.1f',lme_biom(s,1))); hold on;
    text(5,2.4,sprintf('%2.1f',lme_biom(s,2))); hold on;
    text(5,0.6,sprintf('%2.1f',lme_biom(s,3))); hold on;
    text(2,2.2,sprintf('%2.1e',(lme_cF(s,1)+lme_cF(s,2)))); hold on;
    text(4,2.2,sprintf('%2.1e',lme_cP(s,3))); hold on;
    text(3,0.8,sprintf('%2.1e',lme_cD(s,6))); hold on;
    text(4,1.5,sprintf('%2.1e',lme_cD(s,3))); hold on;
    text(5.25,1.5,sprintf('%2.1e',lme_cD(s,4))); hold on;
    xlim([0 6])
    ylim([0 3])
    title(LMEname{s})
    set(gca,'XTickLabel','','YTickLabel','')
    
end
print(f6,'-dpng',[ppath 'Clim_fished_',harv,'_fluxes_Zbio_bent_oligotrophic_LME.png'])


%% Faves of each
for n=1:length(fave)
    s = fave(n);
    
    f7=figure(7);
    subplot(1,3,n)
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',10*lme_prey(s,1)); hold on;
    plot(1,1,'.','color','k','MarkerSize',10*lme_prey(s,4)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',10*lme_biom(s,1)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',10*lme_biom(s,2)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',10*lme_biom(s,3)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',50*(lme_cF(s,1)+lme_cF(s,2))); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',50*lme_cP(s,3)); hold on;
    plot(x3,y1,'color','k','Linewidth',lme_cD(s,6)); hold on;
    if (lme_cD(s,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',10*lme_cD(s,3)); hold on;
    end
    if (lme_cD(s,4) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',10*lme_cD(s,4)); hold on;
    end
    text(1,2.4,sprintf('%2.1f',lme_prey(s,1))); hold on;
    text(1,0.6,sprintf('%2.1f',lme_prey(s,4))); hold on;
    text(3,2.4,sprintf('%2.1f',lme_biom(s,1))); hold on;
    text(5,2.4,sprintf('%2.1f',lme_biom(s,2))); hold on;
    text(5,0.6,sprintf('%2.1f',lme_biom(s,3))); hold on;
    text(2,2.2,sprintf('%2.1e',(lme_cF(s,1)+lme_cF(s,2)))); hold on;
    text(4,2.2,sprintf('%2.1e',lme_cP(s,3))); hold on;
    text(3,0.8,sprintf('%2.1e',lme_cD(s,6))); hold on;
    if (n==2)
        text(4,1.6,sprintf('%2.1f',lme_cD(s,3))); hold on;
        text(5.25,1.5,sprintf('%2.1f',lme_cD(s,4))); hold on;
    else
        text(4,1.6,sprintf('%2.1e',lme_cD(s,3))); hold on;
        text(5.25,1.5,sprintf('%2.1e',lme_cD(s,4))); hold on;
    end
    xlim([0 6])
    ylim([0 3])
    title(LMEname{s})
    set(gca,'XTickLabel','','YTickLabel','')
    
    f8=figure(8);
    subplot(3,1,n)
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',10*lme_prey(s,1)); hold on;
    plot(1,1,'.','color','k','MarkerSize',10*lme_prey(s,4)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',10*lme_biom(s,1)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',10*lme_biom(s,2)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',10*lme_biom(s,3)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',50*(lme_cF(s,1)+lme_cF(s,2))); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',50*lme_cP(s,3)); hold on;
    plot(x3,y1,'color','k','Linewidth',lme_cD(s,6)); hold on;
    if (lme_cD(s,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',10*lme_cD(s,3)); hold on;
    end
    if (lme_cD(s,4) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',10*lme_cD(s,4)); hold on;
    end
    text(1,2.4,sprintf('%2.1f',lme_prey(s,1))); hold on;
    text(1,0.6,sprintf('%2.1f',lme_prey(s,4))); hold on;
    text(3,2.4,sprintf('%2.1f',lme_biom(s,1))); hold on;
    text(5,2.4,sprintf('%2.1f',lme_biom(s,2))); hold on;
    text(5,0.6,sprintf('%2.1f',lme_biom(s,3))); hold on;
    text(2,2.2,sprintf('%2.1e',(lme_cF(s,1)+lme_cF(s,2)))); hold on;
    text(4,2.2,sprintf('%2.1e',lme_cP(s,3))); hold on;
    text(3,0.8,sprintf('%2.1e',lme_cD(s,6))); hold on;
    if (n==2)
        text(4,1.5,sprintf('%2.1f',lme_cD(s,3))); hold on;
        text(5.25,1.5,sprintf('%2.1f',lme_cD(s,4))); hold on;
    else
        text(4,1.5,sprintf('%2.1e',lme_cD(s,3))); hold on;
        text(5.25,1.5,sprintf('%2.1e',lme_cD(s,4))); hold on;
    end
    xlim([0 6])
    ylim([0 3])
    title(LMEname{s})
    set(gca,'XTickLabel','','YTickLabel','')
    
end
print(f7,'-dpng',[ppath 'Clim_fished_',harv,'_fluxes_Zbio_bent_faves_LMEsH.png'])
print(f8,'-dpng',[ppath 'Clim_fished_',harv,'_fluxes_Zbio_bent_faves_LMEsV.png'])

