% Visualize output of POEM Climatology at single locations
% 150 years, monthly means saved

clear all
close all

datap = '/Volumes/FEISTY/NC/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
spots=spots';

load('cmap_ppt_angles.mat')
cmap3=cmap_ppt([5,1,3],:);
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

%% load output
h = 25;         % coeff on Cmax
gam = 4;       % coeff on search area

cfile = ['Dc_Lam700_enc',num2str(gam),'-b200_m400-b175-k086_c',...
    num2str(h),'-b250_D080_A050_nmort1_BE08_noCC_RE00100'];
sfile = ['/Volumes/FEISTY/NC/Matlab_new_size/',cfile,...
    '/Climatol_locs_1meso_All_fish03.mat'];
load(sfile);

fpath = [figp cfile '/'];
if (~isfolder(fpath))
    mkdir(fpath)
end

close all
%%
all_mean=NaN*ones(3,4,length(spots));
z = NaN*ones(length(spots),3);

SP = S_Sml_p;
SF = S_Sml_f;
SD = S_Sml_d;
MP = S_Med_p;
MF = S_Med_f;
MD = S_Med_d;
LP = S_Lrg_p;
LD = S_Lrg_d;
CO = S_Cobalt;

t=1:size(SP,1);
lyr=t((end-12+1):end);

%% Final mean biomass in each size

Fmean=[sf_mean,mf_mean];
Pmean=[sp_mean,mp_mean,lp_mean];
Dmean=[sd_mean,md_mean,ld_mean];
Bmean = b_mean;

all_mean(1:2,1,:) = Fmean';
all_mean(:,2,:) = Pmean';
all_mean(:,3,:) = Dmean';
all_mean(1,4,:) = Bmean';

% Size spectrum (sum stages)
spec = nansum(all_mean,2);

%% Feeding level = con/cmax
SP_lev=squeeze(nanmean(SP(lyr,20,:)))';
SF_lev=squeeze(nanmean(SF(lyr,20,:)))';
SD_lev=squeeze(nanmean(SD(lyr,20,:)))';
MP_lev=squeeze(nanmean(MP(lyr,20,:)))';
MF_lev=squeeze(nanmean(MF(lyr,20,:)))';
MD_lev=squeeze(nanmean(MD(lyr,20,:)))';
LP_lev=squeeze(nanmean(LP(lyr,20,:)))';
LD_lev=squeeze(nanmean(LD(lyr,20,:)))';

Plev=[SP_lev;MP_lev;LP_lev];
Flev=[SF_lev;MF_lev];
Dlev=[SD_lev;MD_lev;LD_lev];

%% Fraction bent biom, zoop losses, zoo biom consumed
z(:,1) = squeeze(nanmean(CO(lyr,1,:) ./ CO(lyr,2,:)));
z(:,2) = squeeze(nanmean(CO(lyr,3,:)));
z(:,3) = squeeze(nanmean(CO(lyr,4,:)));


%% Gross growth efficiency (= nu/consump)
SP_gge=squeeze(nanmean(SP(lyr,15,:)./SP(lyr,14,:)))';
SF_gge=squeeze(nanmean(SF(lyr,15,:)./SF(lyr,14,:)))';
SD_gge=squeeze(nanmean(SD(lyr,15,:)./SD(lyr,14,:)))';
MP_gge=squeeze(nanmean(MP(lyr,15,:)./MP(lyr,14,:)))';
MF_gge=squeeze(nanmean(MF(lyr,15,:)./MF(lyr,14,:)))';
MD_gge=squeeze(nanmean(MD(lyr,15,:)./MD(lyr,14,:)))';
LP_gge=squeeze(nanmean(LP(lyr,15,:)./LP(lyr,14,:)))';
LD_gge=squeeze(nanmean(LD(lyr,15,:)./LD(lyr,14,:)))';

Pgge=[SP_gge;MP_gge;LP_gge];
Fgge=[SF_gge;MF_gge];
Dgge=[SD_gge;MD_gge;LD_gge];

save(sfile,...
    'Pmean','Fmean','Dmean','Bmean','all_mean','z',...
    'Pgge','Fgge','Dgge','Plev','Flev','Dlev','-append');

mlev = [Flev;Plev;Dlev];

%%
for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    
    %% TIME SERIES ----------------------------------------------------
    y=t/12;
    
    F = SF(:,1,s)+MF(:,1,s);
    P = SP(:,1,s)+MP(:,1,s)+LP(:,1,s);
    D = SD(:,1,s)+MD(:,1,s)+LD(:,1,s);
    B = S_Cobalt(:,1,s);
    
    figure(51)
    clf
    plot(y,log10(B),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(F),'r','Linewidth',2); hold on;
    plot(y,log10(P),'b','Linewidth',2); hold on;
    plot(y,log10(D),'k','Linewidth',2); hold on;
    legend('B','F','P','D')
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    ylim([-5 3])
    xlabel('Year')
    ylabel('log10 Biomass (g m^-^2)')
    title(['Clim ' loc])
    stamp(cfile)
    print('-dpng',[fpath 'Ts_1meso_Logmean_biomass_types_' loc '.png'])
    
    %  Comparisons ----------------------------------------------------
    
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
    ylim([-2 3])
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
    
    %% Fraction zoop losses consumed
    f5 = figure(5);
    subplot(4,4,s)
    bar(z(s,:),'k'); hold on;
    xlim([0 4])
    ylim([0 1.1])
    set(gca,'XTick',1:3,'XTickLabel',{'B','Zl','Zb',})
    if (s==5)
        ylabel('Fraction consumed')
    end
    title(loc)
    if (s==3)
        stamp(cfile)
    end
    
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
    ylim([-0.2 0.8])
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
    ylim([-0.2 0.8])
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
print(f21,'-dpng',[fpath 'Clim_All_locs_1meso_Logmean_biomass_axes.png'])
print(f2,'-dpng',[fpath 'Clim_All_locs_1meso_con_level.png'])
print(f5,'-dpng',[fpath 'Clim_All_locs_1meso_frac_zoop_loss.png'])
print(f15,'-dpng',[fpath 'Clim_All_locs_1meso_gge.png'])

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
ylim([-2 3])
set(gca,'XTick',1:16,'XTickLabel',[])
for n=1:16
    text(n,-2.2,spots{n},'HorizontalAlignment','center')
end
ylabel('log10 Mean Biom (g m^-^2) in final year')
title('All stages')
stamp(cfile)
print('-dpng',[fpath 'Clim_All_locs_1meso_tot_mean_biomass_type.png'])

sumspec = squeeze(nansum(nansum(all_mean)));

figure(17);
plot(1:16,log10(sumspec),'k.','MarkerSize',25); hold on;
xlim([0 17])
ylim([-1 3])
set(gca,'XTick',1:16,'XTickLabel',[])
for n=1:16
    text(n,-1.1,spots{n},'HorizontalAlignment','center')
end
ylabel('log10 Mean Biom (g m^-^2) in final year')
title('All fishes and stages')
stamp(cfile)
print('-dpng',[fpath 'Clim_All_locs_1meso_tot_mean_biomass_spec.png'])




