% Visualize output of FEISTY Climatology at single locations
% 150 years, monthly means saved

clear all
close all

warning off 

%% location info
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'MZ','LZ','B','SF','MF','SP','MP','LP','SD','MD','LD'};
cols=cols';
spots=spots';

%% new run with param struct & fish struct
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_A075_Sm025_nmort1_BE08_noCC_RE00100';
fpath = ['/Volumes/FEISTY/NC/Clim_comp_tests/' cfile '/NNU_fishvec_'];
load([fpath,'Climatol_All_fish03_locs.mat'])
figp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/clim_complete/post_proc/pp_figs/',cfile,'/NNU_fishvec_'];

%all outputs time x locs x cols

%% colors
load('cmap_ppt_angles.mat')
cmap3=cmap_ppt([5,1,3],:);
cm=[1 0.5 0;...   %orange
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

set(groot,'defaultAxesColorOrder',cm);

%% param info
mass = param.wc;
mass = repmat(mass,length(spots),1);

A = 4.39;
fc = 0.2;
f0 = 0.6;
epsassim = 0.7;
n = 3/4;
w = logspace(-3, 5);
AvailEnergy = A*w.^n;
Consumption = A / (epsassim*(f0-fc)) * w.^n;
%Andersen & Beyer mortality rate per year (natural + predation)
%physiol mort * growth constant * M^-0.25
AB = (0.35 .* 4.5 .* mass.^(-0.25)) ./365;

stages={'SF','MF','SP','MP','LP','SD','MD','LD'};

%% Take means
all_mean=NaN*ones(3,4,length(spots));
z = NaN*ones(length(spots),3);

[nt,nid,nf] = size(S_bio);
t=1:nt;
lyr=t((end-12+1):end);

%% Final mean biomass in each size

bio_mean=squeeze(mean(S_bio(lyr,:,:)));
prey_mean=squeeze(mean(S_Cobalt(lyr,:,:)));

Fmean=bio_mean(:,param.ix1(1):param.ix2(1));
Pmean=bio_mean(:,param.ix1(2):param.ix2(2));
Dmean=bio_mean(:,param.ix1(3):param.ix2(3));
Bmean = prey_mean(:,1);

all_mean(1:2,1,:) = Fmean';
all_mean(:,2,:) = Pmean';
all_mean(:,3,:) = Dmean';
all_mean(1,4,:) = Bmean';

% Size spectrum (sum stages)
spec = nansum(all_mean,2);

%% Growth rate (nu = energy for biomass production)
mgr_mean=squeeze(mean(S_nu(lyr,:,:)));

Fmgr=mgr_mean(:,param.ix1(1):param.ix2(1))';
Pmgr=mgr_mean(:,param.ix1(2):param.ix2(2))';
Dmgr=mgr_mean(:,param.ix1(3):param.ix2(3))';

%% Consump per biomass (I)
con_mean=squeeze(mean(S_I(lyr,:,:)));

Fcon=con_mean(:,param.ix1(1):param.ix2(1))';
Pcon=con_mean(:,param.ix1(2):param.ix2(2))';
Dcon=con_mean(:,param.ix1(3):param.ix2(3))';

%% Feeding level = con/cmax
lev_mean=squeeze(mean(S_clev(lyr,:,:)));

Flev=lev_mean(:,param.ix1(1):param.ix2(1))';
Plev=lev_mean(:,param.ix1(2):param.ix2(2))';
Dlev=lev_mean(:,param.ix1(3):param.ix2(3))';

%% Fraction zoop losses consumed
z = prey_mean(:,3:5);

%% Production (= nu * biom)
prod_mean=squeeze(mean(S_prod(lyr,:,:)));

Fprod=prod_mean(:,param.ix1(1):param.ix2(1))';
Pprod=prod_mean(:,param.ix1(2):param.ix2(2))';
Dprod=prod_mean(:,param.ix1(3):param.ix2(3))';

%% Reproduction
rep_mean=squeeze(mean(S_rep(lyr,:,:)));

Frep(1,:)=rep_mean(:,param.ix2(1))';
Prep(1,:)=rep_mean(:,param.ix2(2))';
Drep(1,:)=rep_mean(:,param.ix2(3))';
Frep(2,:)=rep_mean(:,param.ix2(1))' .* bio_mean(:,param.ix2(1))';
Prep(2,:)=rep_mean(:,param.ix2(2))' .* bio_mean(:,param.ix2(2))';
Drep(2,:)=rep_mean(:,param.ix2(3))' .* bio_mean(:,param.ix2(3))';

%% Metabolism
met_mean=squeeze(mean(S_met(lyr,:,:)));

Fmet=met_mean(:,param.ix1(1):param.ix2(1))';
Pmet=met_mean(:,param.ix1(2):param.ix2(2))';
Dmet=met_mean(:,param.ix1(3):param.ix2(3))';

%% Predation
pred_mean=squeeze(mean(S_pred(lyr,:,:)));

Fpred=pred_mean(:,param.ix1(1):param.ix2(1))';
Ppred=pred_mean(:,param.ix1(2):param.ix2(2))';
Dpred=pred_mean(:,param.ix1(3):param.ix2(3))';

%% Natural mortality
nat_mean=squeeze(mean(S_nmort(lyr,:,:)));

Fnat=nat_mean(:,param.ix1(1):param.ix2(1))';
Pnat=nat_mean(:,param.ix1(2):param.ix2(2))';
Dnat=nat_mean(:,param.ix1(3):param.ix2(3))';

%% Fishing
fish_mean=squeeze(mean(S_caught(lyr,:,:)));
fish_tot=squeeze(sum(S_caught(lyr,:,:)));
frate_mean=squeeze(mean(S_fmort(lyr,:,:)));

Ffish=fish_mean(:,param.ix1(1):param.ix2(1))';
Pfish=fish_mean(:,param.ix1(2):param.ix2(2))';
Dfish=fish_mean(:,param.ix1(3):param.ix2(3))';

Ffrate=frate_mean(:,param.ix1(1):param.ix2(1))';
Pfrate=frate_mean(:,param.ix1(2):param.ix2(2))';
Dfrate=frate_mean(:,param.ix1(3):param.ix2(3))';

Ftotcatch=fish_tot(:,param.ix1(1):param.ix2(1))';
Ptotcatch=fish_tot(:,param.ix1(2):param.ix2(2))';
Dtotcatch=fish_tot(:,param.ix1(3):param.ix2(3))';

%% Total mortality w/o fishing
Fmort = Fpred + Fnat;
Pmort = Ppred + Pnat;
Dmort = Dpred + Dnat;

%% Total mortality w/ fishing
Fmortf = Fpred + Fnat + Ffish;
Pmortf = Ppred + Pnat + Pfish;
Dmortf = Dpred + Dnat + Dfish;

%% Gross growth efficiency (= nu/consump)
gge = S_nu ./ S_I;
gge_mean=squeeze(mean(gge(lyr,:,:)));

Fgge=gge_mean(:,param.ix1(1):param.ix2(1))';
Pgge=gge_mean(:,param.ix1(2):param.ix2(2))';
Dgge=gge_mean(:,param.ix1(3):param.ix2(3))';

%% Save
save([fpath,'Means_Climatol_All_fish03_locs_lastyr_means.mat'],...
    'Pmean','Fmean','Dmean','all_mean',...
    'Pmgr','Fmgr','Dmgr','Pcon','Fcon','Dcon','z','Pprod','Fprod','Dprod',...
    'Prep','Frep','Drep','Pmet','Fmet','Dmet','Ppred','Fpred','Dpred',...
    'Pnat','Fnat','Dnat','Pfish','Ffish','Dfish','Ptotcatch','Ftotcatch',...
    'Dtotcatch','Pgge','Fgge','Dgge','Plev','Flev','Dlev','Bmean',...
    'Pfrate','Ffrate','Dfrate'); %...
    % 'conF','conP','conD','conZm','conZl','conB'),

%% Plots
mlev = [Flev;Plev;Dlev];
F=squeeze(sum(S_bio(:,:,param.ix1(1):param.ix2(1)),3));
P=squeeze(sum(S_bio(:,:,param.ix1(2):param.ix2(2)),3));
D=squeeze(sum(S_bio(:,:,param.ix1(3):param.ix2(3)),3));

for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    
    %% TIME SERIES ----------------------------------------------------
    y=t/12;
    figure(50)
    clf
    plot(y,log10(squeeze(S_bio(:,s,4:end))),'Linewidth',1); hold on;
    legend('SF','MF','SP','MP','LP','SD','MD','LD')
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    ylim([-5 4])
    xlabel('Year')
    ylabel('log10 Biomass (g m^-^2)')
    title(['Clim ' loc])
    stamp(cfile)
    print('-dpng',[figp 'timeseries_logmean_biomass_' loc '.png'])
    
    figure(51)
    clf
    plot(y,log10(F(:,s)),'r','Linewidth',2); hold on;
    plot(y,log10(P(:,s)),'b','Linewidth',2); hold on;
    plot(y,log10(D(:,s)),'k','Linewidth',2); hold on;
    legend('F','P','D')
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    ylim([-5 5])
    xlabel('Year')
    ylabel('log10 Biomass (g m^-^2)')
    title(['Clim ' loc])
    stamp(cfile)
    print('-dpng',[figp 'timeseries_Logmean_biomass_types_' loc '.png'])
    
    %  TIME SERIES ----------------------------------------------------
    
    %% Biomass
    f21 = figure(21);
    subplot(4,4,s)
    plot(0.5:2:5.5,log10(squeeze(all_mean(:,1,s))),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:2:6,log10(squeeze(all_mean(:,2,s))),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1.5:2:6.5,log10(squeeze(all_mean(:,3,s))),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 6])
    ylim([-3 2])
    set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
    if (s==5)
        ylabel('log10 Mean Biom (g m^-^2) in final year')
    end
    title(loc)
    if (s==3)
        stamp(cfile)
    end
    
    %% Feeding level
    f2=figure(2);
    subplot(4,4,s)
    bar(mlev(:,s),'k')
    ylim([0 1])
    xlim([0 9])
    set(gca,'XTickLabel',[]);
    for n=1:8
        text(n-0.5,-0.2,stages{n},'Rotation',45)
    end
    title(spots{s})
    if (s==5)
        ylabel('Feeding level')
    end
    if (s==16)
        stamp(cfile)
    end
    
    %% Growth rate (nu = energy for biomass production)
    f3 = figure(3);
    subplot(2,2,1)
    plot(s-0.25,Fmgr(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,Pmgr(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,Dmgr(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0,spots{n},'Rotation',45)
    end
    ylabel('Mean growth rate (g g^-^1 d^-^1)')
    title('S')
    
    subplot(2,2,2)
    plot(s-0.25,(Fmgr(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(Pmgr(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dmgr(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0,spots{n},'Rotation',45)
    end
    ylabel('Mean growth/repro rate (g g^-^1 d^-^1)')
    title('M')
    
    subplot(2,2,3)
    plot(s,(Pmgr(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dmgr(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0,spots{n},'Rotation',45)
    end
    ylabel('Mean repro rate (g g^-^1 d^-^1) in final year')
    title('L')
    if (s==3)
        stamp(cfile)
    end
    
    %% Fraction zoop losses consumed
    f5 = figure(5);
    subplot(4,4,s)
    bar(z(s,:),'k'); hold on;
    xlim([0 4])
    ylim([0 1])
    set(gca,'XTick',1:3,'XTickLabel',{'MZ','LZ','Bent'})
    if (s==5)
        ylabel('Fraction flux consumed')
    end
    title(loc)
    if (s==3)
        stamp(cfile)
    end
    
    %% Production (= nu * biom)
    f8 = figure(8);
    subplot(2,2,1)
    plot(s-0.25,Fprod(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,Pprod(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,Dprod(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    ylim([0 0.1])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,-0.01,spots{n},'Rotation',45)
    end
    %ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
    title('S')
    
    subplot(2,2,2)
    plot(s-0.25,(Fprod(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(Pprod(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dprod(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    ylim([0 0.1])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,-0.01,spots{n},'Rotation',45)
    end
    %ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
    title('M')
    
    subplot(2,2,3)
    plot(s,(Pprod(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dprod(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    ylim([0 0.05])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,-0.01,spots{n},'Rotation',45)
    end
    ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
    title('L')
    if (s==1)
        stamp(cfile)
    end
    
    %% Reproduction
    f9 = figure(9);
    subplot(1,2,1)
    plot(s-0.25,Frep(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,Prep(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,Drep(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    ylim([0 0.02])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,-0.001,spots{n},'Rotation',45)
    end
    ylabel('Mean repro rate (g g^-^1 d^-^1) in final year')
    
    subplot(1,2,2)
    plot(s-0.25,(Frep(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(Prep(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Drep(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    ylim([0 0.1])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,-0.005,spots{n},'Rotation',45)
    end
    ylabel('Mean biom reproduced (g d^-^1) in final year')
    if (s==1)
        stamp(cfile)
    end
    
    %% Metabolism
    f10 = figure(10);
    subplot(2,2,1)
    plot(s-0.25,Fmet(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,Pmet(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,Dmet(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0,spots{n},'Rotation',45)
    end
    %ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
    title('S')
    
    subplot(2,2,2)
    plot(s-0.25,(Fmet(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(Pmet(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dmet(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0,spots{n},'Rotation',45)
    end
    %ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
    title('M')
    
    subplot(2,2,3)
    plot(s,(Pmet(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dmet(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0,spots{n},'Rotation',45)
    end
    ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
    title('L')
    if (s==3)
        stamp(cfile)
    end
    
    %% Predation
    f11 = figure(11);
    subplot(1,2,1)
    plot(s-0.25,Fpred(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,Ppred(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,Dpred(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,-0.001,spots{n},'Rotation',45)
    end
    ylabel('Mean predation rate (d^-^1) in final year')
    title('S')
    
    subplot(1,2,2)
    plot(s-0.25,(Fpred(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(Ppred(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dpred(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,-0.001,spots{n},'Rotation',45)
    end
    ylabel('Mean predation rate (d^-^1) in final year')
    title('M')
    if (s==3)
        stamp(cfile)
    end
    
    %% Fishing
    f12 = figure(12);
    subplot(1,2,1)
    plot(s-0.25,Ftotcatch(2,s),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,Ptotcatch(2,s),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,Dtotcatch(2,s),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,-0.1,spots{n},'Rotation',45)
    end
    ylabel('Total catch (g) in final year')
    title('Medium')
    
    subplot(1,2,2)
    plot(s,Ptotcatch(3,s),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,Dtotcatch(3,s),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,-0.1,spots{n},'Rotation',45)
    end
    ylabel('Total catch (g) in final year')
    title('Large')
    if (s==1)
        stamp(cfile)
    end
    
    %% Total mortality w/o fishing
    Fmort = Fpred + Fnat;
    Pmort = Ppred + Pnat;
    Dmort = Dpred + Dnat;
    
    f13=figure(13);
    subplot(2,2,1)
    plot(s-0.25,Fmort(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,Pmort(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,Dmort(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0,spots{n},'Rotation',45)
    end
    if(s==16)
        plot(0:17,AB(1)*ones(18,1),'--k'); hold on;
    end
    ylabel('Mean mortality rate w/o fishing (d^-^1) in final year')
    title('S')
    
    subplot(2,2,2)
    plot(s-0.25,(Fmort(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(Pmort(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dmort(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0,spots{n},'Rotation',45)
    end
    if(s==16)
        plot(0:17,AB(2)*ones(18,1),'--k'); hold on;
    end
    title('M')
    
    subplot(2,2,3)
    plot(s,(Pmort(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dmort(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0,spots{n},'Rotation',45)
    end
    if(s==16)
        plot(0:17,AB(3)*ones(18,1),'--k'); hold on;
    end
    title('L')
    
    %% Total mortality w/ fishing
    Fmortf = Fpred + Fnat + Ffish;
    Pmortf = Ppred + Pnat + Pfish;
    Dmortf = Dpred + Dnat + Dfish;
    
    f14=figure(14);
    subplot(2,2,1)
    plot(s-0.25,Fmortf(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,Pmortf(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,Dmortf(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0,spots{n},'Rotation',45)
    end
    if(s==16)
        plot(0:17,AB(1)*ones(18,1),'--k'); hold on;
    end
    ylabel('Mean mortality rate w/fishing (d^-^1) in final year')
    title('S')
    
    subplot(2,2,2)
    plot(s-0.25,(Fmortf(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(Pmortf(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dmortf(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0,spots{n},'Rotation',45)
    end
    if(s==16)
        plot(0:17,AB(2)*ones(18,1),'--k'); hold on;
    end
    title('M')
    
    subplot(2,2,3)
    plot(s,(Pmortf(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dmortf(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0,spots{n},'Rotation',45)
    end
    if(s==16)
        plot(0:17,AB(3)*ones(18,1),'--k'); hold on;
    end
    title('L')
    
    %% Gross growth efficiency (= nu/consump)
    f15 = figure(15);
    subplot(2,2,1)
    plot(s-0.25,Fgge(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,Pgge(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,Dgge(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    ylim([-0.2 0.8])
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,-0.3,spots{n},'Rotation',45)
    end
    ylabel('Mean gross growth efficiency in final year')
    title('S')
    
    subplot(2,2,2)
    plot(s-0.25,(Fgge(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(Pgge(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dgge(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    ylim([-0.2 0.6])
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,-0.3,spots{n},'Rotation',45)
    end
    %ylabel('Mean gross growth efficiency in final year')
    title('M')
    
    subplot(2,2,3)
    plot(s,(Pgge(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dgge(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    ylim([-0.2 0.5])
    xlim([0 17])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,-0.25,spots{n},'Rotation',45)
    end
    %ylabel('Mean gross growth efficiency in final year')
    title('L')
    if (s==3)
        stamp(cfile)
    end
    
end
print(f21,'-dpng',[figp 'All_oneloc_Logmean_biomass_axes.png'])
print(f2,'-dpng',[figp 'All_oneloc_con_level.png'])
print(f3,'-dpng',[figp 'All_oneloc_nu.png'])
print(f5,'-dpng',[figp 'All_oneloc_frac_zoop_loss.png'])
print(f8,'-dpng',[figp 'All_oneloc_prod.png'])
print(f9,'-dpng',[figp 'All_oneloc_rep.png'])
print(f10,'-dpng',[figp 'All_oneloc_met.png'])
print(f11,'-dpng',[figp 'All_oneloc_pred.png'])
print(f12,'-dpng',[figp 'All_oneloc_catch.png'])
print(f13,'-dpng',[figp 'All_oneloc_mort_nof.png'])
print(f14,'-dpng',[figp 'All_oneloc_mort_f.png'])
print(f15,'-dpng',[figp 'All_oneloc_gge.png'])

%% Sum mean biom over stages
fishsp = squeeze(nansum(all_mean));

figure(16);
plot((1-0.2):16,log10(fishsp(1,:)),'sk','MarkerFaceColor',cmap_ppt(3,:),...
    'MarkerSize',15); hold on;
plot(1:16,log10(fishsp(2,:)),'sk','MarkerFaceColor',cmap_ppt(1,:),...
    'MarkerSize',15); hold on;
plot((1+0.2):17,log10(fishsp(3,:)),'sk','MarkerFaceColor',cmap_ppt(2,:),...
    'MarkerSize',15); hold on;
xlim([0 17])
ylim([-2 4])
set(gca,'XTick',1:16,'XTickLabel',[])
for n=1:16
    text(n,-2.2,spots{n},'HorizontalAlignment','center')
end
ylabel('log10 Mean Biom (g m^-^2) in final year')
title('All stages')
stamp(cfile)
print('-dpng',[figp 'All_oneloc_tot_mean_biomass_type.png'])

sumspec = squeeze(nansum(nansum(all_mean)));

figure(17);
plot(1:16,log10(sumspec),'k.','MarkerSize',25); hold on;
xlim([0 17])
ylim([0 5])
set(gca,'XTick',1:16,'XTickLabel',[])
for n=1:16
    text(n,-2.1,spots{n},'HorizontalAlignment','center')
end
ylabel('log10 Mean Biom (g m^-^2) in final year')
title('All fishes and stages')
stamp(cfile)
print('-dpng',[figp 'All_oneloc_tot_mean_biomass_spec.png'])

%% Fishing rate
figure(54);
plot((1-0.2):16,Ffrate(2,:)*365,'sk','MarkerFaceColor',cmap_ppt(3,:),...
    'MarkerSize',15); hold on;
plot(1:16,Pfrate(3,:)*365,'sk','MarkerFaceColor',cmap_ppt(1,:),...
    'MarkerSize',15); hold on;
plot((1+0.2):17,Dfrate(3,:)*365,'sk','MarkerFaceColor',cmap_ppt(2,:),...
    'MarkerSize',15); hold on;
xlim([0 17])
ylim([0 1])
set(gca,'XTick',1:16,'XTickLabel',[])
for n=1:16
    text(n,-0.05,spots{n},'HorizontalAlignment','center')
end
ylabel('Mean fishing rate (yr^-^1) in final year')
title('Adult stages')
stamp(cfile)
print('-dpng',[figp 'All_oneloc_mean_frate_type.png'])


%% Consump g/g/d --> g/d --> g/y
Fcon = Fcon' .* mass(:,param.ix1(1):param.ix2(1)) .* 365;
Pcon = Pcon' .* mass(:,param.ix1(2):param.ix2(2)) .* 365;
Dcon = Dcon' .* mass(:,param.ix1(3):param.ix2(3)) .* 365;

f18 = figure(18);
subplot(3,1,1)
plot(0.75:1:15.75,Fcon(:,1),'sk',...
    'MarkerFaceColor',cmap_ppt(3,:),...
    'MarkerSize',15); hold on;
plot(1:1:16,Pcon(:,1),'sk',...
    'MarkerFaceColor',cmap_ppt(1,:),...
    'MarkerSize',15); hold on;
plot(1.25:1:16.25,Dcon(:,1),'sk',...
    'MarkerFaceColor',cmap_ppt(2,:),...
    'MarkerSize',15); hold on;
xlim([0 17])
ylabel('S')
set(gca,'XTick',1:16,'XTickLabel',spots)
title('Mean consumption (g y^-^1) in final year')

subplot(3,1,2)
plot(0.75:1:15.75,Fcon(:,2),'sk',...
    'MarkerFaceColor',cmap_ppt(3,:),...
    'MarkerSize',15); hold on;
plot(1:1:16,Pcon(:,2),'sk',...
    'MarkerFaceColor',cmap_ppt(1,:),...
    'MarkerSize',15); hold on;
plot(1.25:1:16.25,Dcon(:,2),'sk',...
    'MarkerFaceColor',cmap_ppt(2,:),...
    'MarkerSize',15); hold on;
xlim([0 17])
ylabel('M')
set(gca,'XTick',1:16,'XTickLabel',spots)

subplot(3,1,3)
plot(1:1:16,Pcon(:,3),'sk',...
    'MarkerFaceColor',cmap_ppt(1,:),...
    'MarkerSize',15); hold on;
plot(1.25:1:16.25,Dcon(:,3),'sk',...
    'MarkerFaceColor',cmap_ppt(2,:),...
    'MarkerSize',15); hold on;
xlim([0 17])
ylabel('L')
set(gca,'XTick',1:16,'XTickLabel',spots)
xlabel('Location')
stamp(cfile)

print(f18,'-dpng',[figp 'All_oneloc_consump_gyr.png'])

%% Consump vs. weight
f19=figure(19);
for s=1:length(spots)
    subplot(2,2,1)
    loglog((mass(s,param.ix1(1):param.ix2(1))),(Fcon(s,:)),'.',...
        'Color',cm(s,:),'MarkerSize',25); hold on;
    title('F')
    xlabel('Mass (g)')
    ylabel('Mean consumption (g y^-^1) in final year')
    %axis([-5 5 -1 5])
    
    subplot(2,2,2)
    loglog((mass(s,param.ix1(2):param.ix2(2))),(Pcon(s,:)),'.',...
        'Color',cm(s,:),'MarkerSize',25); hold on;
    title('P')
    xlabel('Mass (g)')
    %axis([-5 5 -1 5])
    
    subplot(2,2,3)
    loglog((mass(s,param.ix1(3):param.ix2(3))),(Dcon(s,:)),'.',...
        'Color',cm(s,:),'MarkerSize',25); hold on;
    title('D')
    xlabel('Mass (g)')
    %axis([-5 5 -1 5])
    
    subplot(2,2,4)
    loglog((mass(s,param.ix1(3):param.ix2(3))),(Dcon(s,:)),'.',...
        'Color',cm(s,:),'MarkerSize',25); hold on;
    xlabel('Mass (g)')
    legend(spots)
    legend('location','northwest','orientation','horizontal')
    axis([-5 1 -1 1])
    stamp(cfile)
end
subplot(2,2,1)
loglog(w, Consumption,'k')

subplot(2,2,2)
loglog(w, Consumption,'k')

subplot(2,2,3)
loglog(w, Consumption,'k')

subplot(2,2,4)
legend('location','eastoutside')
legend('orientation','vertical')
print(f19,'-dpng',[figp 'All_oneloc_consump_gyr_vs_weight_compare.png'])


