% Visualize output of POEM biol rate eq tests
% Spinup at one location
% 150 years, monthly means saved
% Independently change coeffs for cmax and met fns

clear all
close all

%GFDL/NC/Matlab_new_size/Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100/param_sens/
%cfile = 'Dc_enc50-b210_m4-b175-k060_c50-b250_D075_J075_A050_Sm025_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

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
spots{16} = 'PUP';

nfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,'/param_sens/'];
load([nfile 'Locs_Climatol_All_fish03_means_aenc_acmax_search.mat']);
pname = 'Climatol_All_fish03_means_aenc_acmax_search';


%%
aep = 10:10:100;
acp = 10:10:100;
nj=length(aep);
nk=length(acp);

% Biomass of each type
allF = SF+MF;
allP = SP+MP+LP;
allD = SD+MD+LD;
allB = BI;

% Gut fullness
Fflev = (fSF+fMF)/2;
Pflev = (fSP+fMP+fLP)/3;
Dflev = (fSD+fMD+fLD)/3;
Sflev = (fSF+fSP+fSD)/3;
Mflev = (fMF+fMP+fMD)/3;
Lflev = (fLP+fLD)/2;
Tflev = (fSF+fMF+fSP+fMP+fLP+fSD+fMD+fLD)/7;

% GGE
Fgge = (gSF+gMF)/2;
Pgge = (gSP+gMP+gLP)/3;
Dgge = (gSD+gMD+gLD)/3;
Sgge = (gSF+gSP+gSD)/3;
Mgge = (gMF+gMP+gMD)/3;
Lgge = (gLP+gLD)/2;
Tgge = (gSF+gMF+gSP+gMP+gLP+gSD+gMD+gLD)/7;

%%
FPrat = squeeze(allF./(allF+allP));
DPrat = squeeze(allD./(allD+allP));

jays = [aep 110];
ays = [acp 110];
[agrid,jgrid]=meshgrid(ays,jays);

allF2 = NaN*ones(nj+1,nk+1,16);
allP2 = NaN*ones(nj+1,nk+1,16);
allD2 = NaN*ones(nj+1,nk+1,16);
FP2 = NaN*ones(nj+1,nk+1,16);
DP2 = NaN*ones(nj+1,nk+1,16);
Tgge2 = NaN*ones(nj+1,nk+1,16);
Tflev2 = NaN*ones(nj+1,nk+1,16);
Sgge2 = NaN*ones(nj+1,nk+1,16);
Sflev2 = NaN*ones(nj+1,nk+1,16);
Mgge2 = NaN*ones(nj+1,nk+1,16);
Mflev2 = NaN*ones(nj+1,nk+1,16);
Lgge2 = NaN*ones(nj+1,nk+1,16);
Lflev2 = NaN*ones(nj+1,nk+1,16);

allF2(1:nj,1:nk,:)=allF;
allP2(1:nj,1:nk,:)=allP;
allD2(1:nj,1:nk,:)=allD;
FP2(1:nj,1:nk,:)=FPrat;
DP2(1:nj,1:nk,:)=DPrat;
Tgge2(1:nj,1:nk,:)=Tgge;
Tflev2(1:nj,1:nk,:)=Tflev;
Sgge2(1:nj,1:nk,:)=Sgge;
Sflev2(1:nj,1:nk,:)=Sflev;
Mgge2(1:nj,1:nk,:)=Mgge;
Mflev2(1:nj,1:nk,:)=Mflev;
Lgge2(1:nj,1:nk,:)=Lgge;
Lflev2(1:nj,1:nk,:)=Lflev;

%%
% colors
cmBP=cbrewer('seq','BuPu',50);
cmYOR=cbrewer('seq','YlOrRd',50);

%% Only use 3 domain examples: EBS(10), PUp(16), HOT(13)
sid = [10;16;13];
for s=1:length(sid)
    domain = sid(s);
    loc = spots{domain};
    lname = [loc '_'];
    
    %% Sum mean biom over stages
    f1=figure(1);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(log10(allF2(:,:,domain))))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmBP)
    caxis([-2 2])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'log10 Mean F Biom (g m^-^2)'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(log10(allP2(:,:,domain))))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmBP)
    caxis([-2 2])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('log10 Mean P Biom (g m^-^2)')
    end
    if (s==1)
        ylabel('a_E')
    end
    
    subplot(3,3,s+6)
    pcolor(agrid,jgrid,squeeze(log10(allD2(:,:,domain))))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmBP)
    caxis([-2 2])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('log10 Mean D Biom (g m^-^2)')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    %stamp(pname)
    
    
    f2=figure(2);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(FP2(:,:,domain)))
    cmocean('balance')
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'Fraction F/(F+P)'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(DP2(:,:,domain)))
    colorbar('Position',[0.92 0.5 0.025 0.3],'orientation','vertical')
    cmocean('balance')
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('Fraction D/(D+P)')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    %stamp(pname)
    
    %% Feeding level
    f3=figure(3);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(Tflev2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0.5 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'Mean feeding level'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
    % GGE
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(Tgge2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('Mean gross growth efficiency')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    %stamp(pname)
    
    %% Feeding level only
    f4=figure(4);
    % Small
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(Sflev2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'S mean feeding level'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
    % Medium
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(Mflev2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('M mean feeding level')
    end
    if (s==1)
        ylabel('a_E')
    end
    
    % Large
    subplot(3,3,s+6)
    pcolor(agrid,jgrid,squeeze(Lflev2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('L mean feeding level')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    %stamp(pname)
    
    %% GGE only
    f5=figure(5);
    % Small
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(Sgge2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'S mean gross growth efficiency'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
    % Medium
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(Mgge2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('M mean gross growth efficiency')
    end
    if (s==1)
        ylabel('a_E')
    end
    
    % Large
    subplot(3,3,s+6)
    pcolor(agrid,jgrid,squeeze(Lgge2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('L mean gross growth efficiency')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    %stamp(pname)
    
end %spots
%%
print(f1,'-dpng',[figp pname '_mean_biomass_type_all_3locs.png'])
print(f2,'-dpng',[figp pname '_frac_all_3locs.png'])
print(f3,'-dpng',[figp pname '_flev_gge_all_3locs.png'])
print(f4,'-dpng',[figp pname '_flev_size_all_3locs.png'])
print(f5,'-dpng',[figp pname '_gge_size_all_3locs.png'])



