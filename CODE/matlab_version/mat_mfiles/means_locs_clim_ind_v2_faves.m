% Visualize output of POEM Climatology at single locations
% Foodweb diagram
% 150 years, monthly means saved

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev','T');
sites = T{:,1};
spots = abbrev;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
red = [6,10,12:16];
spots=spots(red)';
shelf = 1:2;
olig = 3:5;
upw = 6:7;
fave = [2;7;4];

dp = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
sname = 'Clim_';
harv = 'All_fish03';
dpath = [datap char(dp) '/'];
fpath = [figp char(dp) '/'];
if (~isdir([figp char(dp)]))
    mkdir([figp char(dp)])
end
cfile = char(dp);
load([dpath sname 'locs_' harv '.mat'])
load([dpath sname 'locs_' harv '_lastyr_sum_mean_biom.mat']);

% Colors
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
%load('/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
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

M_s = 10^((log10(0.001)+log10(0.5))/2);
M_m = 10^((log10(0.5)+log10(250))/2);
M_l = 10^((log10(250)+log10(125000))/2);

%! Body lengths (mm)
% Convert from mm to cm and use their const coeff = 0.01g/cm3
L_s = 10.0 * (M_s/0.01)^(1/3); % small
L_m = 10.0 * (M_m/0.01)^(1/3); % medium
L_l = 10.0 * (M_l/0.01)^(1/3); % large

mass = [M_s;M_m;M_l];
mass = repmat(mass,1,length(spots));
L = [L_s;L_m;L_l];

stages={'SF','MF','SP','MP','LP','SD','MD','LD'};

%% POEM means

mlev = [Flev;Plev;Dlev];
F = squeeze(nansum(all_mean(:,1,:)));
P = squeeze(nansum(all_mean(:,2,:)));
D = squeeze(nansum(all_mean(:,3,:)));
B = squeeze(nansum(all_mean(:,4,:)));
conZ = conZm + conZl;

%% Zoop, det, bent
cpath = ['/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/'];
load([cpath 'cobalt_zoop_biom_means.mat'],'mz_mean_clim','lz_mean_clim','mzloss_mean_clim','lzloss_mean_clim')
load([cpath 'cobalt_det_biom_means.mat'],'det_mean_clim')

gpath='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([gpath 'clim_npp_Dmeans_Ytot.mat'])

load(['/Volumes/GFDL/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%ESM2.6 in mg C m-2 or mg C m-2 d-1
%from mg C m-2 to g(WW) m-2
% 1e-3 g C in 1 mg C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)

z_mean = (mz_mean_clim + lz_mean_clim) * 1e-3 * 9.0;
z_loss = (mzloss_mean_clim+lzloss_mean_clim) * 1e-3 * 9.0;

z_mean_grid = z_mean(ID);
z_loss_grid = z_loss(ID);

det_grid = det_mean_clim(ID) * 1e-3 * 9.0;

mnpp = npp_mean_clim(ID) * 1e-3 * 9.0;

z_mean_locs = z_mean_grid(ids);
z_loss_locs = z_loss_grid(ids);
det_locs = det_grid(ids);
npp_locs = mnpp(ids);

%% Table
Tab=table(sites,npp_locs,z_mean_locs,z_loss_locs,det_locs,B,F,P,D,...
    'VariableNames',{'Location','NPP','Z','Zloss','Det','Bent','F','P','D'});
writetable(Tab,[dpath sname harv '_locs_biomasses.csv'],'Delimiter',',');
save([dpath sname harv '_locs_biomasses.mat'],'Tab');


%% Plots
x0=linspace(-0.8,0.8,9);
x1=linspace(1.2,2.8,9);
x2=linspace(3.2,4.8,9);
x3=linspace(1.2,4.8,9);
x4=5*ones(9,1);
y1=ones(9,1);
y2=2*ones(9,1);
y3=linspace(1.8,1.2,9);


bios(:,1) = npp_locs(red);
bios(:,2) = z_mean_locs(red);
bios(:,3) = B(red);
bios(:,4) = F(red);
bios(:,5) = P(red);
bios(:,6) = D(red);

flux(:,1) = det_locs(red);
flux(:,2) = z_loss_locs(red);
flux(:,3) = conZ(red,1);
flux(:,4) = conF(red,2);
flux(:,5) = conF(red,3);
flux(:,6) = conP(red,3);
flux(:,7) = conB(red,3);

Lbios = log10(bios);
Lflux = log10(flux+eps);

Sbios = (bios - min(bios(:))) / (max(bios(:)) - min(bios(:)));
Sbios2 = 1 + (bios - min(bios(:)))*(200-1) / (max(bios(:)) - min(bios(:)));
SLbios = (Lbios - min(Lbios(:))) / (max(Lbios(:)) - min(Lbios(:)));
Sflux = (flux - min(flux(:))) / (max(flux(:)) - min(flux(:)));
Sflux2 = 1+ (flux - min(flux(:)))*(15-1) / (max(flux(:)) - min(flux(:)));
SLflux = (Lflux - min(Lflux(:))) / (max(Lflux(:)) - min(Lflux(:)));

Sbios = Sbios+eps;
SLbios = SLbios+eps;

for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    n=red(s);
    
    %% Bubble/arrow plot ----------------------------------------------------
    
    figure(1)
    clf
    plot(-1,2,'.','color',[0.5 0 1],'MarkerSize',Sbios2(s,1)); hold on;
    plot(x0,y2,'color',[0.5 0 1],'Linewidth',Sflux2(s,2)); hold on;
    plot(x0,y3,'color',[0.5 0 1],'Linewidth',Sflux2(s,1)); hold on;
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',Sbios2(s,2)); hold on;
    plot(1,1,'.','color','k','MarkerSize',Sbios2(s,3)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',Sbios2(s,4)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',Sflux2(s,3)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',Sflux2(s,4)); hold on;
    plot(x3,y1,'color','k','Linewidth',Sflux2(s,7)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',Sbios2(s,5)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',Sbios2(s,6)); hold on;
    if (conF(n,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',Sflux2(s,5)); hold on;
        text(4,1.5,sprintf('%2.1e',conF(n,3))); hold on;
    else
        text(4,1.5,sprintf('%2.1f',conF(n,3))); hold on;
    end
    if (conP(n,3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',Sflux2(s,6)); hold on;
        text(5.25,1.5,sprintf('%2.1e',conP(n,3))); hold on;
    else
        text(5.25,1.5,sprintf('%2.1f',conP(n,3))); hold on;
    end
    text(-1,2.4,sprintf('%2.1f',npp_locs(n))); hold on;
    text(0,2.2,sprintf('%2.1f',z_loss_locs(n))); hold on;
    text(1,2.4,sprintf('%2.1f',z_mean_locs(n))); hold on;
    if (det_locs(n) >= 0.1)
        text(0.5,1.5,sprintf('%2.1f',det_locs(n))); hold on;
    else
        text(0.5,1.5,sprintf('%2.1e',det_locs(n))); hold on;
    end
    text(1,0.6,sprintf('%2.1f',B(n))); hold on;
    text(3,2.4,sprintf('%2.1f',F(n))); hold on;
    text(5,2.4,sprintf('%2.1f',P(n))); hold on;
    text(5,0.6,sprintf('%2.1f',D(n))); hold on;
    text(2,2.2,sprintf('%2.1e',conZ(n,1))); hold on;
    text(4,2.2,sprintf('%2.1e',conF(n,2))); hold on;
    text(3,0.8,sprintf('%2.1e',conB(n,3))); hold on;
    xlim([-2 6])
    ylim([0 3])
    title(loc)
    set(gca,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','')
    print('-dpng',[fpath sname harv '_fluxes_npp_Zbio_bent_' loc '_v2_reduced.png'])
    
end

%% Subplots of each type
% Shelf Sea
for s=1:length(shelf)
    loc = spots{shelf(s)};
    lname = [loc '_'];
    n=red(shelf(s));
    
    %% Bubble/arrow plot ----------------------------------------------------
    
    f3=figure(3);
    subplot(1,2,s)
    plot(-1,2,'.','color',[0.5 0 1],'MarkerSize',Sbios2(shelf(s),1)); hold on;
    plot(x0,y2,'color',[0.5 0 1],'Linewidth',Sflux2(shelf(s),2)); hold on;
    plot(x0,y3,'color',[0.5 0 1],'Linewidth',Sflux2(shelf(s),1)); hold on;
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',Sbios2(shelf(s),2)); hold on;
    plot(1,1,'.','color','k','MarkerSize',Sbios2(shelf(s),3)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',Sbios2(shelf(s),4)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',Sflux2(shelf(s),3)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',Sflux2(shelf(s),4)); hold on;
    plot(x3,y1,'color','k','Linewidth',Sflux2(shelf(s),7)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',Sbios2(shelf(s),5)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',Sbios2(shelf(s),6)); hold on;
    if (conF(n,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',Sflux2(shelf(s),5)); hold on;
        text(4,1.5,sprintf('%2.1e',conF(n,3))); hold on;
    else
        text(4,1.5,sprintf('%2.1f',conF(n,3))); hold on;
    end
    if (conP(n,3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',Sflux2(shelf(s),6)); hold on;
        text(5.25,1.5,sprintf('%2.1e',conP(n,3))); hold on;
    else
        text(5.25,1.5,sprintf('%2.1f',conP(n,3))); hold on;
    end
    text(-1,2.4,sprintf('%2.1f',npp_locs(n))); hold on;
    text(0,2.2,sprintf('%2.1f',z_loss_locs(n))); hold on;
    text(1,2.4,sprintf('%2.1f',z_mean_locs(n))); hold on;
    if (det_locs(n) >= 0.1)
        text(0.5,1.5,sprintf('%2.1f',det_locs(n))); hold on;
    else
        text(0.5,1.5,sprintf('%2.1e',det_locs(n))); hold on;
    end
    text(1,0.6,sprintf('%2.1f',B(n))); hold on;
    text(3,2.4,sprintf('%2.1f',F(n))); hold on;
    text(5,2.4,sprintf('%2.1f',P(n))); hold on;
    text(5,0.6,sprintf('%2.1f',D(n))); hold on;
    text(2,2.2,sprintf('%2.1e',conZ(n,1))); hold on;
    text(4,2.2,sprintf('%2.1e',conF(n,2))); hold on;
    text(3,0.8,sprintf('%2.1e',conB(n,3))); hold on;
    xlim([-2 6])
    ylim([0 3])
    title(loc)
    set(gca,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','')
    
end
print(f3,'-dpng',[fpath sname harv '_fluxes_npp_Zbio_bent_shelf_reduced.png'])

%% Upwelling
for s=1:length(upw)
    loc = spots{upw(s)};
    lname = [loc '_'];
    n=red(upw(s));
    
    %% Bubble/arrow plot ----------------------------------------------------
    
    f4=figure(4);
    subplot(1,2,s)
    plot(-1,2,'.','color',[0.5 0 1],'MarkerSize',Sbios2(upw(s),1)); hold on;
    plot(x0,y2,'color',[0.5 0 1],'Linewidth',Sflux2(upw(s),2)); hold on;
    plot(x0,y3,'color',[0.5 0 1],'Linewidth',Sflux2(upw(s),1)); hold on;
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',Sbios2(upw(s),2)); hold on;
    plot(1,1,'.','color','k','MarkerSize',Sbios2(upw(s),3)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',Sbios2(upw(s),4)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',Sflux2(upw(s),3)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',Sflux2(upw(s),4)); hold on;
    plot(x3,y1,'color','k','Linewidth',Sflux2(upw(s),7)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',Sbios2(upw(s),5)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',Sbios2(upw(s),6)); hold on;
    if (conF(n,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',Sflux2(upw(s),5)); hold on;
        text(4,1.5,sprintf('%2.1e',conF(n,3))); hold on;
    else
        text(4,1.5,sprintf('%2.1f',conF(n,3))); hold on;
    end
    if (conP(n,3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',Sflux2(upw(s),6)); hold on;
        text(5.25,1.5,sprintf('%2.1e',conP(n,3))); hold on;
    else
        text(5.25,1.5,sprintf('%2.1f',conP(n,3))); hold on;
    end
    text(-1,2.4,sprintf('%2.1f',npp_locs(n))); hold on;
    text(0,2.2,sprintf('%2.1f',z_loss_locs(n))); hold on;
    text(1,2.4,sprintf('%2.1f',z_mean_locs(n))); hold on;
    if (det_locs(n) >= 0.1)
        text(0.5,1.5,sprintf('%2.1f',det_locs(n))); hold on;
    else
        text(0.5,1.5,sprintf('%2.1e',det_locs(n))); hold on;
    end
    text(1,0.6,sprintf('%2.1f',B(n))); hold on;
    text(3,2.4,sprintf('%2.1f',F(n))); hold on;
    text(5,2.4,sprintf('%2.1f',P(n))); hold on;
    text(5,0.6,sprintf('%2.1f',D(n))); hold on;
    text(2,2.2,sprintf('%2.1e',conZ(n,1))); hold on;
    text(4,2.2,sprintf('%2.1e',conF(n,2))); hold on;
    text(3,0.8,sprintf('%2.1e',conB(n,3))); hold on;
    xlim([-2 6])
    ylim([0 3])
    title(loc)
    set(gca,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','')
    
end
print(f4,'-dpng',[fpath sname harv '_fluxes_npp_Zbio_bent_upwelling_reduced.png'])

% Oligotrophic
for s=1:length(olig)
    loc = spots{olig(s)};
    lname = [loc '_'];
    n=red(olig(s));
    
    %% Bubble/arrow plot ----------------------------------------------------
    
    f5=figure(5);
    subplot(2,2,s)
    plot(-1,2,'.','color',[0.5 0 1],'MarkerSize',Sbios2(olig(s),1)); hold on;
    plot(x0,y2,'color',[0.5 0 1],'Linewidth',Sflux2(olig(s),2)); hold on;
    plot(x0,y3,'color',[0.5 0 1],'Linewidth',Sflux2(olig(s),1)); hold on;
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',Sbios2(olig(s),2)); hold on;
    plot(1,1,'.','color','k','MarkerSize',Sbios2(olig(s),3)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',Sbios2(olig(s),4)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',Sflux2(olig(s),3)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',Sflux2(olig(s),4)); hold on;
    plot(x3,y1,'color','k','Linewidth',Sflux2(olig(s),7)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',Sbios2(olig(s),5)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',Sbios2(olig(s),6)); hold on;
    if (conF(n,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',Sflux2(olig(s),5)); hold on;
        text(4,1.5,sprintf('%2.1e',conF(n,3))); hold on;
    else
        text(4,1.5,sprintf('%2.1f',conF(n,3))); hold on;
    end
    if (conP(n,3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',Sflux2(olig(s),6)); hold on;
        text(5.25,1.5,sprintf('%2.1e',conP(n,3))); hold on;
    else
        text(5.25,1.5,sprintf('%2.1f',conP(n,3))); hold on;
    end
    text(-1,2.4,sprintf('%2.1f',npp_locs(n))); hold on;
    text(0,2.2,sprintf('%2.1f',z_loss_locs(n))); hold on;
    text(1,2.4,sprintf('%2.1f',z_mean_locs(n))); hold on;
    if (det_locs(n) >= 0.1)
        text(0.5,1.5,sprintf('%2.1f',det_locs(n))); hold on;
    else
        text(0.5,1.5,sprintf('%2.1e',det_locs(n))); hold on;
    end
    text(1,0.6,sprintf('%2.1f',B(n))); hold on;
    text(3,2.4,sprintf('%2.1f',F(n))); hold on;
    text(5,2.4,sprintf('%2.1f',P(n))); hold on;
    text(5,0.6,sprintf('%2.1f',D(n))); hold on;
    text(2,2.2,sprintf('%2.1e',conZ(n,1))); hold on;
    text(4,2.2,sprintf('%2.1e',conF(n,2))); hold on;
    text(3,0.8,sprintf('%2.1e',conB(n,3))); hold on;
    xlim([-2 6])
    ylim([0 3])
    title(loc)
    set(gca,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','')
    
end
print(f5,'-dpng',[fpath sname harv '_fluxes_npp_Zbio_bent_oligotrophic_reduced.png'])


%% Faves of each
for s=1:length(fave)
    loc = spots{fave(s)};
    lname = [loc '_'];
    n=red(fave(s));
    
    %% Bubble/arrow plot ----------------------------------------------------
    
    f6=figure(6);
    if (s==1)
        subplot('position',[0.04 0.2 0.3 0.6])
    elseif (s==2)
        subplot('position',[0.36 0.2 0.3 0.6])
    else
        subplot('position',[0.68 0.2 0.3 0.6])
    end
    plot(-1,2,'.','color',[0.5 0 1],'MarkerSize',Sbios2(fave(s),1)); hold on;
    plot(x0,y2,'color',[0.5 0 1],'Linewidth',Sflux2(fave(s),2)); hold on;
    plot(x0,y3,'color',[0.5 0 1],'Linewidth',Sflux2(fave(s),1)); hold on;
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',Sbios2(fave(s),2)); hold on;
    plot(1,1,'.','color','k','MarkerSize',Sbios2(fave(s),3)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',Sbios2(fave(s),4)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',Sflux2(fave(s),3)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',Sflux2(fave(s),4)); hold on;
    plot(x3,y1,'color','k','Linewidth',Sflux2(fave(s),7)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',Sbios2(fave(s),5)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',Sbios2(fave(s),6)); hold on;
    if (conF(n,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',Sflux2(fave(s),5)); hold on;
        text(4,1.5,sprintf('%2.1e',conF(n,3)),'HorizontalAlignment','center'); hold on;
    else
        text(4,1.5,sprintf('%2.1f',conF(n,3)),'HorizontalAlignment','center'); hold on;
    end
    if (conP(n,3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',Sflux2(fave(s),6)); hold on;
        text(5.25,1.5,sprintf('%2.1e',conP(n,3)),'HorizontalAlignment','center'); hold on;
    else
        text(5.25,1.5,sprintf('%2.1f',conP(n,3)),'HorizontalAlignment','center'); hold on;
    end
    text(-1,2.6,sprintf('%2.1f',npp_locs(n)),'HorizontalAlignment','center'); hold on;
    text(0,2.3,sprintf('%2.1f',z_loss_locs(n)),'HorizontalAlignment','center'); hold on;
    text(1,2.6,sprintf('%2.1f',z_mean_locs(n)),'HorizontalAlignment','center'); hold on;
    if (det_locs(n) >= 0.1)
        text(0.5,1.5,sprintf('%2.1f',det_locs(n)),'HorizontalAlignment','center'); hold on;
    else
        text(0.5,1.5,sprintf('%2.1e',det_locs(n)),'HorizontalAlignment','center'); hold on;
    end
    text(1,0.6,sprintf('%2.1f',B(n)),'HorizontalAlignment','center'); hold on;
    text(3,2.6,sprintf('%2.1f',F(n)),'HorizontalAlignment','center'); hold on;
    text(5,2.6,sprintf('%2.1f',P(n)),'HorizontalAlignment','center'); hold on;
    text(5,0.6,sprintf('%2.1f',D(n)),'HorizontalAlignment','center'); hold on;
    text(2,2.3,sprintf('%2.1e',conZ(n,1)),'HorizontalAlignment','center'); hold on;
    text(4,2.3,sprintf('%2.1e',conF(n,2)),'HorizontalAlignment','center'); hold on;
    text(3,0.8,sprintf('%2.1e',conB(n,3)),'HorizontalAlignment','center'); hold on;
    xlim([-2 6])
    ylim([0 3])
    title(loc)
    set(gca,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','')
    
    f7=figure(7);
    if (s==1)
        subplot('position',[0.1 0.65 0.8 0.3])
    elseif (s==2)
        subplot('position',[0.1 0.33 0.8 0.3])
    else
        subplot('position',[0.1 0.01 0.8 0.3])
    end
    plot(-1,2,'.','color',[0.5 0 1],'MarkerSize',Sbios2(fave(s),1)); hold on;
    plot(x0,y2,'color',[0.5 0 1],'Linewidth',Sflux2(fave(s),2)); hold on;
    plot(x0,y3,'color',[0.5 0 1],'Linewidth',Sflux2(fave(s),1)); hold on;
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',Sbios2(fave(s),2)); hold on;
    plot(1,1,'.','color','k','MarkerSize',Sbios2(fave(s),3)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',Sbios2(fave(s),4)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',Sflux2(fave(s),3)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',Sflux2(fave(s),4)); hold on;
    plot(x3,y1,'color','k','Linewidth',Sflux2(fave(s),7)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',Sbios2(fave(s),5)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',Sbios2(fave(s),6)); hold on;
    if (conF(n,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',Sflux2(fave(s),5)); hold on;
        text(4,1.5,sprintf('%2.1e',conF(n,3)),'HorizontalAlignment','center'); hold on;
    else
        text(4,1.5,sprintf('%2.1f',conF(n,3)),'HorizontalAlignment','center'); hold on;
    end
    if (conP(n,3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',Sflux2(fave(s),6)); hold on;
        text(5.25,1.5,sprintf('%2.1e',conP(n,3)),'HorizontalAlignment','center'); hold on;
    else
        text(5.25,1.5,sprintf('%2.1f',conP(n,3)),'HorizontalAlignment','center'); hold on;
    end
    text(-1,2.9,sprintf('%2.1f',npp_locs(n)),'HorizontalAlignment','center'); hold on;
    text(0,2.3,sprintf('%2.1f',z_loss_locs(n)),'HorizontalAlignment','center'); hold on;
    text(1,2.9,sprintf('%2.1f',z_mean_locs(n)),'HorizontalAlignment','center'); hold on;
    if (det_locs(n) >= 0.1)
        text(0.5,1.5,sprintf('%2.1f',det_locs(n)),'HorizontalAlignment','center'); hold on;
    else
        text(0.5,1.5,sprintf('%2.1e',det_locs(n)),'HorizontalAlignment','center'); hold on;
    end
    text(1,0.6,sprintf('%2.1f',B(n)),'HorizontalAlignment','center'); hold on;
    text(3,2.9,sprintf('%2.1f',F(n)),'HorizontalAlignment','center'); hold on;
    text(5,2.9,sprintf('%2.1f',P(n)),'HorizontalAlignment','center'); hold on;
    text(5,0.6,sprintf('%2.1f',D(n)),'HorizontalAlignment','center'); hold on;
    text(2,2.3,sprintf('%2.1e',conZ(n,1)),'HorizontalAlignment','center'); hold on;
    text(4,2.3,sprintf('%2.1e',conF(n,2)),'HorizontalAlignment','center'); hold on;
    text(3,0.8,sprintf('%2.1e',conB(n,3)),'HorizontalAlignment','center'); hold on;
    xlim([-2 6])
    ylim([0 3.5])
    text(-1.75,1.5,loc,'FontWeight','bold')
    set(gca,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','')
    
end
print(f6,'-dpng',[fpath sname harv '_fluxes_npp_Zbio_bent_favesH_reduced.png'])
print(f7,'-dpng',[fpath sname harv '_fluxes_npp_Zbio_bent_favesV_reduced.png'])

%% No fluxes
for s=1:length(fave)
    loc = spots{fave(s)};
    lname = [loc '_'];
    n=red(fave(s));
    
    %% Bubble/arrow plot ----------------------------------------------------
    
    f8=figure(8);
    if (s==1)
        subplot('position',[0.04 0.2 0.3 0.6])
    elseif (s==2)
        subplot('position',[0.36 0.2 0.3 0.6])
    else
        subplot('position',[0.68 0.2 0.3 0.6])
    end
    plot(-1,2,'.','color',[0.5 0 1],'MarkerSize',Sbios2(fave(s),1)); hold on;
    plot(x0,y2,'color',[0.5 0 1],'Linewidth',Sflux2(fave(s),2)); hold on;
    plot(x0,y3,'color',[0.5 0 1],'Linewidth',Sflux2(fave(s),1)); hold on;
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',Sbios2(fave(s),2)); hold on;
    plot(1,1,'.','color','k','MarkerSize',Sbios2(fave(s),3)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',Sbios2(fave(s),4)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',Sflux2(fave(s),3)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',Sflux2(fave(s),4)); hold on;
    plot(x3,y1,'color','k','Linewidth',Sflux2(fave(s),7)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',Sbios2(fave(s),5)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',Sbios2(fave(s),6)); hold on;
    if (conF(n,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',Sflux2(fave(s),5)); hold on;
    end
    if (conP(n,3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',Sflux2(fave(s),6)); hold on;
    end
    text(-1,2.4,sprintf('%2.1f',npp_locs(n)),'HorizontalAlignment','center'); hold on;
    text(1,2.4,sprintf('%2.1f',z_mean_locs(n)),'HorizontalAlignment','center'); hold on;
    text(1,0.8,sprintf('%2.1f',B(n)),'HorizontalAlignment','center'); hold on;
    text(3,2.4,sprintf('%2.1f',F(n)),'HorizontalAlignment','center'); hold on;
    text(5,2.4,sprintf('%2.1f',P(n)),'HorizontalAlignment','center'); hold on;
    text(5,0.8,sprintf('%2.1f',D(n)),'HorizontalAlignment','center'); hold on;
    xlim([-2 6])
    ylim([0.5 2.75])
    title(loc)
    set(gca,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','')
    
    f9=figure(9);
    if (s==1)
        subplot('position',[0.1 0.65 0.8 0.3])
    elseif (s==2)
        subplot('position',[0.1 0.33 0.8 0.3])
    else
        subplot('position',[0.1 0.01 0.8 0.3])
    end
    plot(-1,2,'.','color',[0.5 0 1],'MarkerSize',Sbios2(fave(s),1)); hold on;
    plot(x0,y2,'color',[0.5 0 1],'Linewidth',Sflux2(fave(s),2)); hold on;
    plot(x0,y3,'color',[0.5 0 1],'Linewidth',Sflux2(fave(s),1)); hold on;
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',Sbios2(fave(s),2)); hold on;
    plot(1,1,'.','color','k','MarkerSize',Sbios2(fave(s),3)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',Sbios2(fave(s),4)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',Sflux2(fave(s),3)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',Sflux2(fave(s),4)); hold on;
    plot(x3,y1,'color','k','Linewidth',Sflux2(fave(s),7)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',Sbios2(fave(s),5)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',Sbios2(fave(s),6)); hold on;
    if (conF(n,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',Sflux2(fave(s),5)); hold on;
    end
    if (conP(n,3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',Sflux2(fave(s),6)); hold on;
    end
    text(-1,3.0,sprintf('%2.1f',npp_locs(n)),'HorizontalAlignment','center'); hold on;
    text(1,3.0,sprintf('%2.1f',z_mean_locs(n)),'HorizontalAlignment','center'); hold on;
    text(1,0.6,sprintf('%2.1f',B(n)),'HorizontalAlignment','center'); hold on;
    text(3,3.0,sprintf('%2.1f',F(n)),'HorizontalAlignment','center'); hold on;
    text(5,3.0,sprintf('%2.1f',P(n)),'HorizontalAlignment','center'); hold on;
    text(5,0.6,sprintf('%2.1f',D(n)),'HorizontalAlignment','center'); hold on;
    xlim([-2 6])
    ylim([0 3.5])
    text(-1.75,1.5,loc,'FontWeight','bold')
    set(gca,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','')
    
end
print(f8,'-dpng',[fpath sname harv '_fluxes_npp_Zbio_bent_favesH_red_noFlux.png'])
print(f9,'-dpng',[fpath sname harv '_fluxes_npp_Zbio_bent_favesV_red_noFlux.png'])


%% No numbers
for s=1:length(fave)
    loc = spots{fave(s)};
    lname = [loc '_'];
    n=red(fave(s));
    
    %% Bubble/arrow plot ----------------------------------------------------
    
    f10=figure(10);
    if (s==1)
        subplot('position',[0.04 0.2 0.3 0.6])
    elseif (s==2)
        subplot('position',[0.36 0.2 0.3 0.6])
    else
        subplot('position',[0.68 0.2 0.3 0.6])
    end
    plot(-1,2,'.','color',[0.5 0 1],'MarkerSize',Sbios2(fave(s),1)); hold on;
    plot(x0,y2,'color',[0.5 0 1],'Linewidth',Sflux2(fave(s),2)); hold on;
    plot(x0,y3,'color',[0.5 0 1],'Linewidth',Sflux2(fave(s),1)); hold on;
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',Sbios2(fave(s),2)); hold on;
    plot(1,1,'.','color','k','MarkerSize',Sbios2(fave(s),3)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',Sbios2(fave(s),4)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',Sflux2(fave(s),3)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',Sflux2(fave(s),4)); hold on;
    plot(x3,y1,'color','k','Linewidth',Sflux2(fave(s),7)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',Sbios2(fave(s),5)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',Sbios2(fave(s),6)); hold on;
    if (conF(n,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',Sflux2(fave(s),5)); hold on;
    end
    if (conP(n,3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',Sflux2(fave(s),6)); hold on;
    end
    xlim([-2 6])
    ylim([0.75 2.25])
    title(loc)
    set(gca,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','')
    
    f11=figure(11);
    if (s==1)
        subplot('position',[0.1 0.65 0.8 0.3])
    elseif (s==2)
        subplot('position',[0.1 0.33 0.8 0.3])
    else
        subplot('position',[0.1 0.01 0.8 0.3])
    end
    plot(-1,2,'.','color',[0.5 0 1],'MarkerSize',Sbios2(fave(s),1)); hold on;
    plot(x0,y2,'color',[0.5 0 1],'Linewidth',Sflux2(fave(s),2)); hold on;
    plot(x0,y3,'color',[0.5 0 1],'Linewidth',Sflux2(fave(s),1)); hold on;
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',Sbios2(fave(s),2)); hold on;
    plot(1,1,'.','color','k','MarkerSize',Sbios2(fave(s),3)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',Sbios2(fave(s),4)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',Sflux2(fave(s),3)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',Sflux2(fave(s),4)); hold on;
    plot(x3,y1,'color','k','Linewidth',Sflux2(fave(s),7)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',Sbios2(fave(s),5)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',Sbios2(fave(s),6)); hold on;
    if (conF(n,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',Sflux2(fave(s),5)); hold on;
    end
    if (conP(n,3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',Sflux2(fave(s),6)); hold on;
    end
    xlim([-2 6])
    ylim([0.75 2.5])
    text(-1.75,1.5,loc,'FontWeight','bold')
    set(gca,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','')
    
end
print(f10,'-dpng',[fpath sname harv '_fluxes_npp_Zbio_bent_favesH_red_noNums.png'])
print(f11,'-dpng',[fpath sname harv '_fluxes_npp_Zbio_bent_favesV_red_noNums.png'])
