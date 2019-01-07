% Visualize output of POEM biol rate eq tests
% Spinup at one location
% 150 years, monthly means saved
% Independently change coeffs for cmax and met fns

clear all
close all

%GFDL/NC/Matlab_new_size/Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100/param_sens/
%cfile = 'Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A050_Sm025_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_enc70-b200_m4-b175-k063_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
%cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

pp = ['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];
figp = [pp cfile '/param_sens/'];
if (~isdir(figp))
    mkdir(figp)
end

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
spots=spots';
spots{16}='PUP';
%Only use 3 domain examples: EBS(10), PUp(16), HOT(13)

nfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,'/param_sens/'];
load([nfile 'Locs_Climatol_All_fish03_means_bm_bc_search.mat']);
pname = 'Climatol_All_fish03_means_bm_bc_search';


%%
bcp = 0.1:0.025:0.325;
bmp = 0.1:0.025:0.325;
nj=length(bmp);
nk=length(bcp);

% Biomass of each type
allF = SF+MF;
allP = SP+MP+LP;
allD = SD+MD+LD;
allB = BI;

%%
FPrat = squeeze(allF./(allF+allP));
DPrat = squeeze(allD./(allD+allP));

jays = [bmp 0.35];
ays = [bcp 0.35];
[agrid,jgrid]=meshgrid(ays,jays);

allF2 = NaN*ones(nj+1,nk+1,16);
allP2 = NaN*ones(nj+1,nk+1,16);
allD2 = NaN*ones(nj+1,nk+1,16);
FP2 = NaN*ones(nj+1,nk+1,16);
DP2 = NaN*ones(nj+1,nk+1,16);

allF2(1:nj,1:nk,:)=allF;
allP2(1:nj,1:nk,:)=allP;
allD2(1:nj,1:nk,:)=allD;
FP2(1:nj,1:nk,:)=FPrat;
DP2(1:nj,1:nk,:)=DPrat;

%%
% colors
cmBP=cbrewer('seq','BuPu',50);

%% Only use 3 domain examples: EBS(10), PUp(16), HOT(13)
sid = [10;16;13];
for s=1:length(sid)
    domain = sid(s);
    loc = spots{domain};
    lname = [loc '_'];
    
    %% Sum mean biom over stages
    f2=figure(2);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(log10(allF2(:,:,domain))))
    colormap(cmBP)
    caxis([-2 2])
    set(gca,'XTick',bcp(1:2:end)+0.0125,'XTickLabel',{'.1','.15','.2','.25','.3'},...
        'YTick',bmp(1:2:end)+0.0125,'YTickLabel',bmp(1:2:end))
    if (s==2)
        title({loc; 'log10 Mean F Biom (g m^-^2)'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('b_M')
    end
    
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(log10(allP2(:,:,domain))))
    colormap(cmBP)
    caxis([-2 2])
    set(gca,'XTick',bcp(1:2:end)+0.0125,'XTickLabel',{'.1','.15','.2','.25','.3'},...
        'YTick',bmp(1:2:end)+0.0125,'YTickLabel',bmp(1:2:end))
    if (s==2)
        title('log10 Mean P Biom (g m^-^2)')
    end
    if (s==1)
        ylabel('b_M')
    end
    
    subplot(3,3,s+6)
    pcolor(agrid,jgrid,squeeze(log10(allD2(:,:,domain))))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmBP)
    caxis([-2 2])
    set(gca,'XTick',bcp(1:2:end)+0.0125,'XTickLabel',{'.1','.15','.2','.25','.3'},...
        'YTick',bmp(1:2:end)+0.0125,'YTickLabel',bmp(1:2:end))
    if (s==2)
        title('log10 Mean D Biom (g m^-^2)')
    end
    xlabel('b_C')
    if (s==1)
        ylabel('b_M')
    end
    %stamp(pname)
    
    f5=figure(5);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(FP2(:,:,domain)))
    cmocean('balance')
    caxis([0 1])
    set(gca,'XTick',bcp(1:2:end)+0.0125,'XTickLabel',{'.1','.15','.2','.25','.3'},...
        'YTick',bmp(1:2:end)+0.0125,'YTickLabel',bmp(1:2:end))
    if (s==2)
        title({loc; 'Fraction F/(F+P)'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('b_M')
    end
    %stamp(pname)
    
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(DP2(:,:,domain)))
    colorbar('Position',[0.92 0.55 0.025 0.3],'orientation','vertical')
    cmocean('balance')
    caxis([0 1])
    set(gca,'XTick',bcp(1:2:end)+0.0125,'XTickLabel',{'.1','.15','.2','.25','.3'},...
        'YTick',bmp(1:2:end)+0.0125,'YTickLabel',bmp(1:2:end))
    if (s==2)
        title('Fraction D/(D+P)')
    end
    xlabel('b_C')
    if (s==1)
        ylabel('b_M')
    end
    %stamp(pname)
    
end %spots
%%
print(f2,'-dpng',[figp pname '_mean_biomass_type_all_3locs.png'])
print(f5,'-dpng',[figp pname '_frac_all_3locs.png'])

