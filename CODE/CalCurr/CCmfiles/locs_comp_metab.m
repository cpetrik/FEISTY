% Visualize output of POEM Climatology at single locations
% 150 years, monthly means saved

clear all
close all

warning off 

datap = '/Volumes/GFDL/CSV/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
%figp = '/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
spots=spots';

sname = 'Clim_';
harv = 'All_fish03';
fpath = [figp 'Bio_rates/'];

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
%load('/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
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

%% Load
dpA = 'Dc_enc70-b200_m4-b175-k08_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
dpathA = [datap char(dpA) '/'];
load([dpathA sname 'locs_' harv '_lastyr_sum_mean_biom.mat'],...
    'Pmet','Fmet','Dmet');
AFmet = Fmet;
APmet = Pmet;
ADmet = Dmet;
clear Pmet Fmet Dmet

dpB = 'Dc_enc70-b200_m4-b175-k09_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
dpathB = [datap char(dpB) '/'];
load([dpathB sname 'locs_' harv '_lastyr_sum_mean_biom.mat'],...
    'Pmet','Fmet','Dmet');
BFmet = Fmet;
BPmet = Pmet;
BDmet = Dmet;
clear Pmet Fmet Dmet

cfile = 'Dc_enc70-b200_m4-b175_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100_kt08-09';

%%
aDPmet = ADmet./APmet;
bDPmet = BDmet./BPmet;

%%
for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    
    %% Metabolism
    %Medium
    f1 = figure(1);
    subplot(1,2,1)
    plot(s-0.25,(AFmet(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(APmet(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(ADmet(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    ylim([0.002 0.026])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0.001,spots{n},'Rotation',45)
    end
    ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
    title('M kt=0.08')
    
    subplot(1,2,2)
    plot(s-0.25,(BFmet(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(BPmet(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(BDmet(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    ylim([0.002 0.026])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0.001,spots{n},'Rotation',45)
    end
    ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
    title('M kt=0.09')
    if (s==3)
        stamp(cfile)
    end
    
    
    %Large
    f2 = figure(2);
    subplot(1,2,1)
    plot(s,(APmet(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(ADmet(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    ylim([0.5e-3 9e-3])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0,spots{n},'Rotation',45)
    end
    ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
    title('L kt=0.08')
    
    subplot(1,2,2)
    plot(s,(BPmet(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(BDmet(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 17])
    ylim([0.5e-3 9e-3])
    set(gca,'XTickLabel',[]);
    for n=1:16
        text(n-0.5,0,spots{n},'Rotation',45)
    end
    ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
    title('L kt=0.09')
    if (s==3)
        stamp(cfile)
    end
    
end
print(f1,'-dpng',[fpath sname harv '_' cfile '_locs_met_M.png'])
print(f2,'-dpng',[fpath sname harv '_' cfile '_locs_met_L.png'])

%%
for s=1:11
    loc = spots{s};
    lname = [loc '_'];
    
    %% Metabolism
    %Medium
    f3 = figure(3);
    subplot(1,2,1)
    plot(s-0.25,(AFmet(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(APmet(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(ADmet(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 12])
    ylim([0.002 0.008])
    set(gca,'XTickLabel',[]);
    for n=1:11
        text(n-0.5,0.0019,spots{n},'Rotation',45)
    end
    ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
    title('M kt=0.08')
    
    subplot(1,2,2)
    plot(s-0.25,(BFmet(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(BPmet(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(BDmet(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 12])
    ylim([0.002 0.008])
    set(gca,'XTickLabel',[]);
    for n=1:11
        text(n-0.5,0.0019,spots{n},'Rotation',45)
    end
    ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
    title('M kt=0.09')
    if (s==3)
        stamp(cfile)
    end
    
    
    %Large
    f4 = figure(4);
    subplot(1,2,1)
    plot(s,(APmet(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(ADmet(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 12])
    ylim([0.5e-3 2.5e-3])
    set(gca,'XTickLabel',[]);
    for n=1:11
        text(n-0.5,0.4e-3,spots{n},'Rotation',45)
    end
    ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
    title('L kt=0.08')
    
    subplot(1,2,2)
    plot(s,(BPmet(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(BDmet(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 12])
    ylim([0.5e-3 2.5e-3])
    set(gca,'XTickLabel',[]);
    for n=1:11
        text(n-0.5,0.4e-3,spots{n},'Rotation',45)
    end
    ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
    title('L kt=0.09')
    if (s==3)
        stamp(cfile)
    end
    
end
print(f3,'-dpng',[fpath sname harv '_' cfile '_locs_met_M_cool.png'])
print(f4,'-dpng',[fpath sname harv '_' cfile '_locs_met_L_cool.png'])


