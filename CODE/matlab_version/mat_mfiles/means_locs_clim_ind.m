% Visualize output of POEM Climatology at single locations
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
spots=spots';

dp = 'Dc_enc70-b200_cm25_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
sname = 'Clim_';
harv = 'All_fish03_Juve00';
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
npp_locs = npp_grid(ids);

%% Table
Tab=table(sites,npp_locs,z_mean_locs,z_loss_locs,det_locs,B,F,P,D,...
    'VariableNames',{'Location','NPP','Z','Zloss','Det','Bent','F','P','D'});
writetable(Tab,[dpath sname harv '_locs_biomasses.csv'],'Delimiter',',');
save([dpath sname harv '_locs_biomasses.mat'],'Tab');


%% Plots
x1=linspace(1.2,2.8,9);
x2=linspace(3.2,4.8,9);
x3=linspace(1.2,4.8,9);
x4=5*ones(9,1);
y1=ones(9,1);
y2=2*ones(9,1);
y3=linspace(1.8,1.2,9);

for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    
    %% Bubble/arrow plot ----------------------------------------------------
    
    figure(1)
    clf
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',10*z_mean_locs(s)); hold on;
    plot(1,1,'.','color','k','MarkerSize',10*det_locs(s)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',10*F(s)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',10*P(s)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',10*D(s)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',50*conZ(s,1)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',50*conF(s,2)); hold on;
    plot(x3,y1,'color','k','Linewidth',50*conB(s,3)); hold on;
    if (conF(s,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',10*conF(s,3)); hold on;
    end
    if (conP(s,3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',10*conP(s,3)); hold on;
    end
    text(1,2.4,sprintf('%2.1f',z_mean_locs(s))); hold on;
    text(1,0.6,sprintf('%2.1f',det_locs(s))); hold on;
    text(3,2.4,sprintf('%2.1f',F(s))); hold on;
    text(5,2.4,sprintf('%2.1f',P(s))); hold on;
    text(5,0.6,sprintf('%2.1f',D(s))); hold on;
    text(2,2.2,sprintf('%2.1e',conZ(s,1))); hold on;
    text(4,2.2,sprintf('%2.1e',conF(s,2))); hold on;
    text(3,0.8,sprintf('%2.1e',conB(s,3))); hold on;
    text(4,1.5,sprintf('%2.1e',conF(s,3))); hold on;
    text(5.25,1.5,sprintf('%2.1e',conP(s,3))); hold on;
    xlim([0 6])
    ylim([0 3])
    title(loc)
    set(gca,'XTickLabel','','YTickLabel','')
    print('-dpng',[fpath sname harv '_fluxes_Zbio_' loc '.png'])
    
    figure(2)
    clf
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',10*z_loss_locs(s)); hold on;
    plot(1,1,'.','color','k','MarkerSize',10*det_locs(s)); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',10*F(s)); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',10*P(s)); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',10*D(s)); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',50*conZ(s,1)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',50*conF(s,2)); hold on;
    plot(x3,y1,'color','k','Linewidth',50*conB(s,3)); hold on;
    if (conF(s,3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',10*conF(s,3)); hold on;
    end
    if (conP(s,3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',10*conP(s,3)); hold on;
    end
    text(1,2.4,sprintf('%2.1f',z_mean_locs(s))); hold on;
    text(1,0.6,sprintf('%2.1f',det_locs(s))); hold on;
    text(3,2.4,sprintf('%2.1f',F(s))); hold on;
    text(5,2.4,sprintf('%2.1f',P(s))); hold on;
    text(5,0.6,sprintf('%2.1f',D(s))); hold on;
    text(2,2.2,sprintf('%2.1e',conZ(s,1))); hold on;
    text(4,2.2,sprintf('%2.1e',conF(s,2))); hold on;
    text(3,0.8,sprintf('%2.1e',conB(s,3))); hold on;
    text(4,1.5,sprintf('%2.1e',conF(s,3))); hold on;
    text(5.25,1.5,sprintf('%2.1e',conP(s,3))); hold on;
    xlim([0 6])
    ylim([0 3])
    title(loc)
    set(gca,'XTickLabel','','YTickLabel','')
    print('-dpng',[fpath sname harv '_fluxes_Zloss_' loc '.png'])
    
end

%% Subplots of each type
shelf = [6;10];
upw = [15:16];
olig = [12:14];
fave = [6;16;14];

% Shelf Sea
for s=1:length(shelf)
    loc = spots{shelf(s)};
    lname = [loc '_'];
    
    %% Bubble/arrow plot ----------------------------------------------------
    
    f3=figure(3);
    subplot(1,2,s)
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',10*z_mean_locs(shelf(s))); hold on;
    plot(1,1,'.','color','k','MarkerSize',10*det_locs(shelf(s))); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',10*F(shelf(s))); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',10*P(shelf(s))); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',10*D(shelf(s))); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',50*conZ(shelf(s),1)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',50*conF(shelf(s),2)); hold on;
    plot(x3,y1,'color','k','Linewidth',50*conB(shelf(s),3)); hold on;
    if (conF(shelf(s),3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',50*conF(shelf(s),3)); hold on;
    end
    if (conP(shelf(s),3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',50*conP(shelf(s),3)); hold on;
    end
    text(1,2.4,sprintf('%2.1f',z_mean_locs(shelf(s)))); hold on;
    text(1,0.6,sprintf('%2.1f',det_locs(shelf(s)))); hold on;
    text(3,2.4,sprintf('%2.1f',F(shelf(s)))); hold on;
    text(5,2.4,sprintf('%2.1f',P(shelf(s)))); hold on;
    text(5,0.6,sprintf('%2.1f',D(shelf(s)))); hold on;
    text(2,2.2,sprintf('%2.1e',conZ(shelf(s),1))); hold on;
    text(4,2.2,sprintf('%2.1e',conF(shelf(s),2))); hold on;
    text(3,0.8,sprintf('%2.1e',conB(shelf(s),3))); hold on;
    text(4,1.5,sprintf('%2.1e',conF(shelf(s),3))); hold on;
    text(5.25,1.5,sprintf('%2.1e',conP(shelf(s),3))); hold on;
    xlim([0 6])
    ylim([0 3])
    title(loc)
    set(gca,'XTickLabel','','YTickLabel','')
    
end
print(f3,'-dpng',[fpath sname harv '_fluxes_Zbio_shelf.png'])

% Upwelling
for s=1:length(upw)
    loc = spots{upw(s)};
    lname = [loc '_'];
    
    %% Bubble/arrow plot ----------------------------------------------------
    
    f4=figure(4);
    subplot(1,2,s)
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',10*z_mean_locs(upw(s))); hold on;
    plot(1,1,'.','color','k','MarkerSize',10*det_locs(upw(s))); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',10*F(upw(s))); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',10*P(upw(s))); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',10*D(upw(s))); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',50*conZ(upw(s),1)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',50*conF(upw(s),2)); hold on;
    plot(x3,y1,'color','k','Linewidth',50*conB(upw(s),3)); hold on;
    if (conF(upw(s),3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',50*conF(upw(s),3)); hold on;
    end
    if (conP(upw(s),3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',50*conP(upw(s),3)); hold on;
    end
    text(1,2.4,sprintf('%2.1f',z_mean_locs(upw(s)))); hold on;
    text(1,0.6,sprintf('%2.1f',det_locs(upw(s)))); hold on;
    text(3,2.4,sprintf('%2.1f',F(upw(s)))); hold on;
    text(5,2.4,sprintf('%2.1f',P(upw(s)))); hold on;
    text(5,0.6,sprintf('%2.1f',D(upw(s)))); hold on;
    text(2,2.2,sprintf('%2.1e',conZ(upw(s),1))); hold on;
    text(4,2.2,sprintf('%2.1e',conF(upw(s),2))); hold on;
    text(3,0.8,sprintf('%2.1e',conB(upw(s),3))); hold on;
    text(4,1.5,sprintf('%2.1f',conF(upw(s),3))); hold on;
    text(5.25,1.5,sprintf('%2.1f',conP(upw(s),3))); hold on;
    xlim([0 6])
    ylim([0 3])
    title(loc)
    set(gca,'XTickLabel','','YTickLabel','')
    
end
print(f4,'-dpng',[fpath sname harv '_fluxes_Zbio_upwelling.png'])

% Oligotrophic
for s=1:length(olig)
    loc = spots{olig(s)};
    lname = [loc '_'];
    
    %% Bubble/arrow plot ----------------------------------------------------
    
    f5=figure(5);
    subplot(2,2,s)
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',10*z_mean_locs(olig(s))); hold on;
    plot(1,1,'.','color','k','MarkerSize',10*det_locs(olig(s))); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',10*F(olig(s))); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',10*P(olig(s))); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',10*D(olig(s))); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',50*conZ(olig(s),1)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',50*conF(olig(s),2)); hold on;
    plot(x3,y1,'color','k','Linewidth',50*conB(olig(s),3)); hold on;
    if (conF(olig(s),3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',50*conF(olig(s),3)); hold on;
    end
    if (conP(olig(s),3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',50*conP(olig(s),3)); hold on;
    end
    text(1,2.4,sprintf('%2.1f',z_mean_locs(olig(s)))); hold on;
    text(1,0.6,sprintf('%2.1f',det_locs(olig(s)))); hold on;
    text(3,2.4,sprintf('%2.1f',F(olig(s)))); hold on;
    text(5,2.4,sprintf('%2.1f',P(olig(s)))); hold on;
    text(5,0.6,sprintf('%2.1f',D(olig(s)))); hold on;
    text(2,2.2,sprintf('%2.1e',conZ(olig(s),1))); hold on;
    text(4,2.2,sprintf('%2.1e',conF(olig(s),2))); hold on;
    text(3,0.8,sprintf('%2.1e',conB(olig(s),3))); hold on;
    text(4,1.5,sprintf('%2.1f',conF(olig(s),3))); hold on;
    text(5.25,1.5,sprintf('%2.1f',conP(olig(s),3))); hold on;
    xlim([0 6])
    ylim([0 3])
    title(loc)
    set(gca,'XTickLabel','','YTickLabel','')
    
end
print(f5,'-dpng',[fpath sname harv '_fluxes_Zbio_oligotrophic.png'])


%% Faves of each
for s=1:length(fave)
    loc = spots{fave(s)};
    lname = [loc '_'];
    
    %% Bubble/arrow plot ----------------------------------------------------
    
    f6=figure(6);
    subplot(1,3,s)
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',10*z_mean_locs(fave(s))); hold on;
    plot(1,1,'.','color','k','MarkerSize',10*det_locs(fave(s))); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',10*F(fave(s))); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',10*P(fave(s))); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',10*D(fave(s))); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',50*conZ(fave(s),1)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',50*conF(fave(s),2)); hold on;
    plot(x3,y1,'color','k','Linewidth',50*conB(fave(s),3)); hold on;
    if (conF(fave(s),3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',50*conF(fave(s),3)); hold on;
    end
    if (conP(fave(s),3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',50*conP(fave(s),3)); hold on;
    end
    text(1,2.4,sprintf('%2.1f',z_mean_locs(fave(s)))); hold on;
    text(1,0.6,sprintf('%2.1f',det_locs(fave(s)))); hold on;
    text(3,2.4,sprintf('%2.1f',F(fave(s)))); hold on;
    text(5,2.4,sprintf('%2.1f',P(fave(s)))); hold on;
    text(5,0.6,sprintf('%2.1f',D(fave(s)))); hold on;
    text(2,2.2,sprintf('%2.1e',conZ(fave(s),1))); hold on;
    text(4,2.2,sprintf('%2.1e',conF(fave(s),2))); hold on;
    text(3,0.8,sprintf('%2.1e',conB(fave(s),3))); hold on;
    if (s==1)
        text(4,1.6,sprintf('%2.1e',conF(shelf(s),3))); hold on;
        text(5.25,1.5,sprintf('%2.1e',conP(shelf(s),3))); hold on;
    else
        text(4,1.6,sprintf('%2.1f',conF(fave(s),3))); hold on;
        text(5.25,1.5,sprintf('%2.1f',conP(fave(s),3))); hold on;
    end
    xlim([0 6])
    ylim([0 3])
    title(loc)
    set(gca,'XTickLabel','','YTickLabel','')
    
    f7=figure(7);
    subplot(3,1,s)
    plot(1,2,'.','color',[0.5 0.5 0.5],'MarkerSize',10*z_mean_locs(fave(s))); hold on;
    plot(1,1,'.','color','k','MarkerSize',10*det_locs(fave(s))); hold on;
    plot(3,2,'.','color',cmap3(1,:),'MarkerSize',10*F(fave(s))); hold on;
    plot(5,2,'.','color',cmap3(2,:),'MarkerSize',10*P(fave(s))); hold on;
    plot(5,1,'.','color',cmap3(3,:),'MarkerSize',10*D(fave(s))); hold on;
    plot(x1,y2,'color',[0.5 0.5 0.5],'Linewidth',50*conZ(fave(s),1)); hold on;
    plot(x2,y2,'color',cmap3(1,:),'Linewidth',50*conF(fave(s),2)); hold on;
    plot(x3,y1,'color','k','Linewidth',50*conB(fave(s),3)); hold on;
    if (conF(fave(s),3) > 0)
        plot(x2,y3,'color',cmap3(1,:),'Linewidth',50*conF(fave(s),3)); hold on;
    end
    if (conP(fave(s),3) > 0)
        plot(x4,y3,'color',cmap3(2,:),'Linewidth',50*conP(fave(s),3)); hold on;
    end
    text(1,2.4,sprintf('%2.1f',z_mean_locs(fave(s)))); hold on;
    text(1,0.6,sprintf('%2.1f',det_locs(fave(s)))); hold on;
    text(3,2.4,sprintf('%2.1f',F(fave(s)))); hold on;
    text(5,2.4,sprintf('%2.1f',P(fave(s)))); hold on;
    text(5,0.6,sprintf('%2.1f',D(fave(s)))); hold on;
    text(2,2.2,sprintf('%2.1e',conZ(fave(s),1))); hold on;
    text(4,2.2,sprintf('%2.1e',conF(fave(s),2))); hold on;
    text(3,0.8,sprintf('%2.1e',conB(fave(s),3))); hold on;
    if (s==1)
        text(4,1.5,sprintf('%2.1e',conF(shelf(s),3))); hold on;
        text(5.25,1.5,sprintf('%2.1e',conP(shelf(s),3))); hold on;
    else
        text(4,1.5,sprintf('%2.1f',conF(fave(s),3))); hold on;
        text(5.25,1.5,sprintf('%2.1f',conP(fave(s),3))); hold on;
    end
    xlim([0 6])
    ylim([0 3])
    title(loc)
    set(gca,'XTickLabel','','YTickLabel','')
    
end
print(f6,'-dpng',[fpath sname harv '_fluxes_Zbio_favesH.png'])
print(f7,'-dpng',[fpath sname harv '_fluxes_Zbio_favesV.png'])

