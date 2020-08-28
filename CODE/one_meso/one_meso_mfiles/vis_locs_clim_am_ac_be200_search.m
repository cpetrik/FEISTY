% Visualize output of FEISTY biol rate eq tests
% Climatology at one location
% 150 years, monthly means saved
% Independently change coeffs for cmax and enc fns

clear all
close all

pfile = 'Dc_Lam700_enc70-b200_m4-b175-k086_c20-b250_D080_A050_nmort1_BE08_noCC_RE00100';
nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/',pfile,'/param_sens/'];
load([nfile 'Locs_Climatol_1meso_All_fish03_noCC_means_aenc_acmax_be200_search.mat'])

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
figp = [pp pfile '/param_sens/'];
if (~isfolder(figp))
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

pname = 'Climatol_1meso_All_fish03_means_aenc_acmax_be200_search';

%%
aep = 1:3:25;
acp = 9:2:25;
nj=length(aep);
nk=length(acp);

% Biomass of each type
allF = SF+MF;
allP = SP+MP+LP;
allD = SD+MD+LD;
allB = BI;

% Gut fullness
% Fflev = (fSF+fMF)/2;
% Pflev = (fSP+fMP+fLP)/3;
% Dflev = (fSD+fMD+fLD)/3;
% Sflev = (fSF+fSP+fSD)/3;
% Mflev = (fMF+fMP+fMD)/3;
% Lflev = (fLP+fLD)/2;
% Tflev = (fSF+fMF+fSP+fMP+fLP+fSD+fMD+fLD)/7;
%just adults
Fflev = (fMF);
Pflev = (fLP);
Dflev = (fLD);
Tflev = (fMF+fLP+fLD)/3;

% GGE
% Fgge = (gSF+gMF)/2;
% Pgge = (gSP+gMP+gLP)/3;
% Dgge = (gSD+gMD+gLD)/3;
% Sgge = (gSF+gSP+gSD)/3;
% Mgge = (gMF+gMP+gMD)/3;
% Lgge = (gLP+gLD)/2;
% Tgge = (gSF+gMF+gSP+gMP+gLP+gSD+gMD+gLD)/7;
%just adults
Fgge = (gMF);
Pgge = (gLP);
Dgge = (gLD);
Tgge = (gMF+gLP+gLD)/3;

% ZL & ZB

%%
FPrat = squeeze(allF./(allF+allP));
DPrat = squeeze(allD./(allD+allP));

jays = [aep 28];
ays = [acp 27];
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
ZL2 = NaN*ones(nj+1,nk+1,16);
ZB2 = NaN*ones(nj+1,nk+1,16);

allF2(1:nj,1:nk,:)=allF;
allP2(1:nj,1:nk,:)=allP;
allD2(1:nj,1:nk,:)=allD;
FP2(1:nj,1:nk,:)=FPrat;
DP2(1:nj,1:nk,:)=DPrat;

Tgge2(1:nj,1:nk,:)=Tgge;
Tflev2(1:nj,1:nk,:)=Tflev;

Sgge2(1:nj,1:nk,:)=Fgge;
Sflev2(1:nj,1:nk,:)=Fflev;

Mgge2(1:nj,1:nk,:)=Pgge;
Mflev2(1:nj,1:nk,:)=Pflev;

Lgge2(1:nj,1:nk,:)=Dgge;
Lflev2(1:nj,1:nk,:)=Dflev;

ZL2(1:nj,1:nk,:)=ZL;
ZB2(1:nj,1:nk,:)=ZB;

%%
% colors
cmBP=cbrewer('seq','BuPu',10,'PCHIP');
cmYOR=cbrewer('seq','YlOrRd',10,'PCHIP');

%% Only use 3 Pac domain examples: EBS(10), PUp(16), HOT(13)
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
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
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
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
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
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('log10 Mean D Biom (g m^-^2)')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    stamp('')
    
    
    f2=figure(2);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(FP2(:,:,domain)))
    cmocean('balance')
    caxis([0 1])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
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
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('Fraction D/(D+P)')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    stamp('')
    
    %% Feeding level
    f3=figure(3);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(Tflev2(:,:,domain)))
    colorbar('Position',[0.92 0.675 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
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
    colorbar('Position',[0.92 0.325 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([-0.5 0.5])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('Mean gross growth efficiency')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    stamp('')
    
    %% Feeding level only
    f4=figure(4);
    % Small
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(Sflev2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title({loc; 'MF mean feeding level'})
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
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('LP mean feeding level')
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
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('LD mean feeding level')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    stamp('')
    
    %% GGE only
    f5=figure(5);
    % Small
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(Sgge2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    cmocean('balance')
    caxis([-0.5 0.5])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title({loc; 'MF mean gross growth efficiency'})
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
    cmocean('balance')
    caxis([-0.5 0.5])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('LP mean gross growth efficiency')
    end
    if (s==1)
        ylabel('a_E')
    end
    
    % Large
    subplot(3,3,s+6)
    pcolor(agrid,jgrid,squeeze(Lgge2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    cmocean('balance')
    caxis([-0.5 0.5])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('LD mean gross growth efficiency')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    stamp('')
    
    %% Zooplankton consumed
    f6=figure(6);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(ZL2(:,:,domain)))
    colormap(cmYOR)
    %caxis([0.1 1.1])
    caxis([0.5 5.5])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title({loc; 'Frac ZL con'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(ZB2(:,:,domain)))
    colorbar('Position',[0.92 0.5 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0.5 5.5])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('Frac ZB con')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    stamp('')
    
end %spots
%%
print(f1,'-dpng',[figp pname '_mean_biomass_type_all_3locsPac.png'])
print(f2,'-dpng',[figp pname '_frac_all_3locsPac.png'])
print(f3,'-dpng',[figp pname '_flev_gge_all_3locsPac.png'])
print(f4,'-dpng',[figp pname '_flev_type_all_3locsPac.png'])
print(f5,'-dpng',[figp pname '_gge_type_all_3locsPac.png'])
print(f6,'-dpng',[figp pname '_zcon_all_3locsPac.png'])


%% Use 3 non-Pac domain examples: NS(6), EEP(15), BATS(14)
sid = [6;15;14];

for s=1:length(sid)
    domain = sid(s);
    loc = spots{domain};
    lname = [loc '_'];
    
    %% Sum mean biom over stages
    f11=figure(11);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(log10(allF2(:,:,domain))))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmBP)
    caxis([-2 2])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
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
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
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
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('log10 Mean D Biom (g m^-^2)')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    stamp('')
    
    
    f12=figure(12);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(FP2(:,:,domain)))
    cmocean('balance')
    caxis([0 1])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
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
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('Fraction D/(D+P)')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    stamp('')
    
    %% Feeding level
    f13=figure(13);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(Tflev2(:,:,domain)))
    colorbar('Position',[0.92 0.675 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
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
    colorbar('Position',[0.92 0.325 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([-0.5 0.5])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('Mean gross growth efficiency')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    stamp('')
    
    %% Feeding level only
    f14=figure(14);
    % Small
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(Sflev2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title({loc; 'MF mean feeding level'})
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
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('LP mean feeding level')
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
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('LD mean feeding level')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    stamp('')
    
    %% GGE only
    f15=figure(15);
    % Small
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(Sgge2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    cmocean('balance')
    caxis([-0.5 0.5])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title({loc; 'MF mean gross growth efficiency'})
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
    cmocean('balance')
    caxis([-0.5 0.5])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('LP mean gross growth efficiency')
    end
    if (s==1)
        ylabel('a_E')
    end
    
    % Large
    subplot(3,3,s+6)
    pcolor(agrid,jgrid,squeeze(Lgge2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    cmocean('balance')
    caxis([-0.5 0.5])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('LD mean gross growth efficiency')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    stamp('')
    
    %% Zooplankton consumed
    f16=figure(16);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(ZL2(:,:,domain)))
    colormap(cmYOR)
    caxis([0.5 5.5])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title({loc; 'Frac ZL con'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(ZB2(:,:,domain)))
    colorbar('Position',[0.92 0.5 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0.5 5.5])
    set(gca,'XTick',(9:3:25)+0.5,'XTickLabel',9:3:25,...
        'YTick',(5:5:25)+5,'YTickLabel',5:5:25)
    if (s==2)
        title('Frac ZB con')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    stamp('')
    
end %spots
%%
print(f11,'-dpng',[figp pname '_mean_biomass_type_all_3locsPac.png'])
print(f12,'-dpng',[figp pname '_frac_all_3locs.png'])
print(f13,'-dpng',[figp pname '_flev_gge_all_3locs.png'])
print(f14,'-dpng',[figp pname '_flev_type_all_3locs.png'])
print(f15,'-dpng',[figp pname '_gge_type_all_3locs.png'])
print(f16,'-dpng',[figp pname '_zcon_all_3locs.png'])




