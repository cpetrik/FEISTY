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

sname = 'Clim_locs_All_fish03';
harv = 'All_fish03';

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

kays = 0:0.25:1;
np = length(kays);
fishsp  = NaN*ones(4,length(spots),np);
mgr     = NaN*ones(3,length(spots),np);
mprod   = NaN*ones(3,length(spots),np);
mrep    = NaN*ones(3,length(spots),np);
mrepb   = NaN*ones(3,length(spots),np);
mgge    = NaN*ones(3,length(spots),np);
for n = 1:length(kays)
    Ka = kays(n);
    tk = num2str(100+int64(10*Ka));
    close all
    dp = ['Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_',...
        'J100_A050_Sm025_nmort1_BE05_noCC_RE00100_K',tk(2:end)];
    dpath = [datap char(dp) '/'];
    fpath = [figp char(dp) '/'];
    cfile = char(dp);
    cfile2 = ['Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_',...
        'A050_Sm025_nmort1_BE05_noCC_RE00100_Jtests'];
    
    load([dpath sname '_lastyr_sum_mean_biom.mat']);
    
    %% Biomass of each type
    fishsp(:,:,n) = squeeze(nansum(all_mean));
    %% Adult Growth rate (nu - energy for biomass production)
    mgr(1,:,n) = Fmgr(2,:);
    mgr(2,:,n) = Pmgr(3,:);
    mgr(3,:,n) = Dmgr(3,:);
    %% Adult Production (= nu * biom)
    mprod(1,:,n) = Fprod(2,:);
    mprod(2,:,n) = Pprod(3,:);
    mprod(3,:,n) = Dprod(3,:);
    %% Reproduction rate
    mrep(1,:,n) = Frep(1,:);
    mrep(2,:,n) = Prep(1,:);
    mrep(3,:,n) = Drep(1,:);
    %% Reproduction biomass
    mrepb(1,:,n) = Frep(2,:);
    mrepb(2,:,n) = Prep(2,:);
    mrepb(3,:,n) = Drep(2,:);
    %% Gross growth efficiency (= nu/consump)
    mgge(1,:,n) = Fgge(1,:);
    mgge(2,:,n) = Pgge(1,:);
    mgge(3,:,n) = Dgge(1,:);
    
end
cfile2 = ['Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_',...
        'J100_A050_Sm025_nmort1_BE05_noCC_RE00100_Ktests'];
% Save values for all locs and all bees that combo
save([datap 'Bio_rates/' cfile2 '.mat'],'fishsp','mgr','mprod','mrep',...
    'mrepb','mgge');


%%
for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    
    %% Sum mean biom over stages
    f1=figure(1);
    subplot(4,4,s)
    plot(1-0.25:np,log10(squeeze(fishsp(1,s,:))),'sk','MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:np,log10(squeeze(fishsp(2,s,:))),'sk','MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1+0.25:np+1,log10(squeeze(fishsp(3,s,:))),'sk','MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 np+1])
    ylim([-2 2])
    if (s==5)
        ylabel('log10 Mean Biom (g m^-^2) in final year')
    end
    set(gca,'XTick',1:np,'XTickLabel',[]);
    for t=1:np
        text(t,-2.1,num2str(kays(t)),'Rotation',45,'HorizontalAlignment','right')
    end
    stamp(cfile2)
    title([loc ' All stages'])
    
    %% Adult Growth rate
    f2=figure(2);
    subplot(4,4,s)
    plot(1-0.25:np,(squeeze(mgr(1,s,:))),'sk','MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:np,(squeeze(mgr(2,s,:))),'sk','MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1+0.25:np+1,(squeeze(mgr(3,s,:))),'sk','MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 np+1])
    ylim([-0.01 0.02])
    if (s==5)
        ylabel('Mean growth rate (g g^-^1 d^-^1)')
    end
    set(gca,'XTick',1:np,'XTickLabel',[]);
    for t=1:np
        text(t,-2.1,num2str(kays(t)),'Rotation',45,'HorizontalAlignment','right')
    end
    stamp(cfile2)
    title([loc ' Adults'])
    
    %% Adult Production
    f3=figure(3);
    subplot(4,4,s)
    plot(1-0.25:np,(squeeze(mprod(1,s,:))),'sk','MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:np,(squeeze(mprod(2,s,:))),'sk','MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1+0.25:np+1,(squeeze(mprod(3,s,:))),'sk','MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 np+1])
    %ylim([-2 2])
    if (s==5)
        ylabel('Mean production rate (g g^-^1 d^-^1)')
    end
    set(gca,'XTick',1:np,'XTickLabel',[]);
    for t=1:np
        text(t,-2.1,num2str(kays(t)),'Rotation',45,'HorizontalAlignment','right')
    end
    stamp(cfile2)
    title([loc ' Adults'])
    
    %% Reproduction rate
    f4=figure(4);
    subplot(4,4,s)
    plot(1-0.25:np,(squeeze(mrep(1,s,:))),'sk','MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:np,(squeeze(mrep(2,s,:))),'sk','MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1+0.25:np+1,(squeeze(mrep(3,s,:))),'sk','MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 np+1])
    ylim([0 0.015])
    if (s==5)
        ylabel('Mean reproduction rate (g g^-^1 d^-^1)')
    end
    set(gca,'XTick',1:np,'XTickLabel',[]);
    for t=1:np
        text(t,-2.1,num2str(kays(t)),'Rotation',45,'HorizontalAlignment','right')
    end
    stamp(cfile2)
    title([loc ' Adults'])
    
    %% Reproduction biomass
    f5=figure(5);
    subplot(4,4,s)
    plot(1-0.25:np,(squeeze(mrepb(1,s,:))),'sk','MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:np,(squeeze(mrepb(2,s,:))),'sk','MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1+0.25:np+1,(squeeze(mrepb(3,s,:))),'sk','MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 np+1])
    %ylim([-2 2])
    if (s==5)
        ylabel('Mean reproduction biom (g d^-^1)')
    end
    set(gca,'XTick',1:np,'XTickLabel',[]);
    for t=1:np
        text(t,-2.1,num2str(kays(t)),'Rotation',45,'HorizontalAlignment','right')
    end
    stamp(cfile2)
    title([loc ' Adults'])
    
    %% Gross growth efficiency
    f6=figure(6);
    subplot(4,4,s)
    plot(1-0.25:np,(squeeze(mgge(1,s,:))),'sk','MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:np,(squeeze(mgge(2,s,:))),'sk','MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1+0.25:np+1,(squeeze(mgge(3,s,:))),'sk','MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 np+1])
    ylim([0.2 0.6])
    if (s==5)
        ylabel('Mean gross growth efficiency (g g^-^1 d^-^1)')
    end
    set(gca,'XTick',1:np,'XTickLabel',[]);
    for t=1:np
        text(t,-2.1,num2str(kays(t)),'Rotation',45,'HorizontalAlignment','right')
    end
    stamp(cfile2)
    title([loc ' Adults'])
    
end %spots
print(f1,'-dpng',[figp sname cfile2 '_tot_mean_biomass_type_all_locs.png'])
print(f2,'-dpng',[figp sname cfile2 '_growth_all_locs.png'])
print(f3,'-dpng',[figp sname cfile2 '_prod_all_locs.png'])
print(f4,'-dpng',[figp sname cfile2 '_rep_all_locs.png'])
print(f5,'-dpng',[figp sname cfile2 '_repbiom_all_locs.png'])
print(f6,'-dpng',[figp sname cfile2 '_gge_all_locs.png'])


