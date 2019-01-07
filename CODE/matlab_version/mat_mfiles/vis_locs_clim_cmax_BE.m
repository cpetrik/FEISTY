% Visualize output of POEM biol rate eq tests
% Climatol at one location
% 150 years, monthly means saved

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

kays = 20:5:30;
bees = 0.05:0.025:0.1;

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
spots=spots';

sname = 'Clim_locs_All_fish03_';
harv = 'All_fish03';

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
load('MyColormaps.mat')
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

AvailEnergy = A*mass.^n;
Consumption = A / (epsassim*(f0-fc)) * mass.^n;

stages={'SF','MF','SP','MP','LP','SD','MD','LD'};

%%
nc = length(kays);
ne = length(bees);
allF  = NaN*ones(length(spots),nc+1,ne+1);
allP  = NaN*ones(length(spots),nc+1,ne+1);
allD  = NaN*ones(length(spots),nc+1,ne+1);
mcon  = NaN*ones(3,nc+1,ne+1);
mlev  = NaN*ones(3,nc+1,ne+1);
for C = 1:length(kays)
    h = kays(C);
    tcfn = num2str(h);
    for E = 1:length(bees)
        beff = bees(E);
        tbe = num2str(100+int64(100*beff));
        
        close all
        dp = ['Dc_enc70-b200_cm',tcfn,'_m-b175-k09_fcrit20_c-b250',...
            '_D075_J100_A050_Sm025_nmort1_BE',tbe(2:end),'_noCC_RE00100'];
        dpath = [datap char(dp) '/'];
        fpath = [figp char(dp) '/'];
        %         if (~isdir([figp char(dp)]))
        %             mkdir([figp char(dp)])
        %         end
        cfile = char(dp);
        cfile2 = ['Dc_enc70_m-b175-k09_fcrit20_D075_J100_A050_Sm025_nmort1_noCC_RE00100_cmaxBEtests'];
        
        load([dpath sname 'lastyr_sum_mean_biom.mat']);
        
        %% Biomass of each type
        allF(:,C,E) = squeeze(nansum(all_mean(:,1,:)));
        allP(:,C,E) = squeeze(nansum(all_mean(:,2,:)));
        allD(:,C,E) = squeeze(nansum(all_mean(:,3,:)));
        
        %% Consump vs. weight
        Pcon = Pcon .* mass .* 365;
        Fcon = Fcon .* mass(1:2,:) .* 365;
        Dcon = Dcon .* mass .* 365;
        Fcon(3,:) = nan;
        Tcon = [Fcon,Pcon,Dcon];
        Mcon = nanmean(Tcon,2);
        mcon(:,C,E) = Mcon;
        
        %% Feeding level
        Flev(3,:) = nan;
        Tlev = [Flev,Plev,Dlev];
        Mlev = nanmean(Tlev,2);
        mlev(:,C,E) = Mlev;
        
    end %BE
end %Cmax

%% Save values for all locs and all bees that combo
save([datap 'Bio_rates/' cfile2 '.mat'],'allF','allP','allD','mcon','mlev');

%%
kays2 = [kays 35];
bees2 = [bees 0.125];
[bgrid,kgrid]=meshgrid(bees2,kays2);
FPrat = allF./(allF+allP);
DPrat = allD./(allD+allP);
for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    
    %% Sum mean biom over stages
    f2=figure(2);
    subplot(4,4,s)
    pcolor(bgrid,kgrid,squeeze(log10(allF(s,:,:))))
    colorbar
    caxis([-2 2])
    set(gca,'XTick',bees,'XTickLabel',bees,...
        'YTick',kays,'YTickLabel',kays)
    if (s==2)
        title({'log10 Mean F Biom (g m^-^2) in final year'; loc})
    else
        title(loc)
    end
    xlabel('Bent eff')
    ylabel('Cmax wgt exp')
    stamp(cfile2)
    
    f3=figure(3);
    subplot(4,4,s)
    pcolor(bgrid,kgrid,squeeze(log10(allP(s,:,:))))
    colorbar
    caxis([-2 2])
    set(gca,'XTick',bees,'XTickLabel',bees,...
        'YTick',kays,'YTickLabel',kays)
    if (s==2)
        title({'log10 Mean P Biom (g m^-^2) in final year'; loc})
    else
        title(loc)
    end
    xlabel('Bent eff')
    ylabel('Cmax wgt exp')
    stamp(cfile2)
    
    
    f4=figure(4);
    subplot(4,4,s)
    pcolor(bgrid,kgrid,squeeze(log10(allD(s,:,:))))
    colorbar
    caxis([-2 2])
    set(gca,'XTick',bees,'XTickLabel',bees,...
        'YTick',kays,'YTickLabel',kays)
    if (s==2)
        title({'log10 Mean D Biom (g m^-^2) in final year'; loc})
    else
        title(loc)
    end
    xlabel('Bent eff')
    ylabel('Cmax wgt exp')
    stamp(cfile2)
    
    f5=figure(5);
    subplot(4,4,s)
    pcolor(bgrid,kgrid,squeeze(FPrat(s,:,:)))
    colorbar
    colormap(cmap_color_rb)
    caxis([0 1])
    set(gca,'XTick',bees,'XTickLabel',bees,...
        'YTick',kays,'YTickLabel',kays)
    if (s==2)
        title({'Frac F:P in final year'; loc})
    else
        title(loc)
    end
    xlabel('Bent eff')
    ylabel('Cmax wgt exp')
    stamp(cfile2)
    
    f6=figure(6);
    subplot(4,4,s)
    pcolor(bgrid,kgrid,squeeze(DPrat(s,:,:)))
    colorbar
    colormap(cmap_color_rb)
    caxis([0 1])
    set(gca,'XTick',bees,'XTickLabel',bees,...
        'YTick',kays,'YTickLabel',kays)
    if (s==2)
        title({'Frac D:P in final year'; loc})
    else
        title(loc)
    end
    xlabel('Bent eff')
    ylabel('Cmax wgt exp')
    stamp(cfile2)
    
end %spots
print(f2,'-dpng',[figp sname cfile2 '_totF_mean_biomass_type_all_locs.png'])
print(f3,'-dpng',[figp sname cfile2 '_totP_mean_biomass_type_all_locs.png'])
print(f4,'-dpng',[figp sname cfile2 '_totD_mean_biomass_type_all_locs.png'])
print(f5,'-dpng',[figp sname cfile2 '_FP_frac_all_locs.png'])
print(f6,'-dpng',[figp sname cfile2 '_DP_frac_all_locs.png'])

%% Feeding level
figure(10)
subplot(2,2,1)
pcolor(bgrid,kgrid,squeeze(mlev(1,:,:)))
colorbar
caxis([0 1])
set(gca,'XTick',bees,'XTickLabel',bees,...
        'YTick',kays,'YTickLabel',kays)
xlabel('Bent eff')
    ylabel('Cmax wgt exp')
title('S Mean feeding level')

subplot(2,2,2)
pcolor(bgrid,kgrid,squeeze(mlev(2,:,:)))
colorbar
caxis([0 1])
set(gca,'XTick',bees,'XTickLabel',bees,...
        'YTick',kays,'YTickLabel',kays)
xlabel('Bent eff')
    ylabel('Cmax wgt exp')
title('M Mean feeding level')

subplot(2,2,3)
pcolor(bgrid,kgrid,squeeze(mlev(3,:,:)))
colorbar
caxis([0 1])
set(gca,'XTick',bees,'XTickLabel',bees,...
        'YTick',kays,'YTickLabel',kays)
xlabel('Bent eff')
    ylabel('Cmax wgt exp')
title('L Mean feeding level')
print('-dpng',[figp sname cfile2 '_mean_flev_size_all_locs.png'])
