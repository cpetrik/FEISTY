% Visualize output of POEM biol rate eq tests
% Spinup at one location
% 150 years, monthly means saved
% Independently change coeffs for cmax and met fns

clear all
close all

%GFDL/NC/Matlab_new_size/Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100/param_sens/
cfile = 'Dc_enc50-b210_m4-b175-k060_c50-b250_D075_J075_A050_Sm025_nmort1_BE08_noCC_RE00100';

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

allF2(1:nj,1:nk,:)=allF;
allP2(1:nj,1:nk,:)=allP;
allD2(1:nj,1:nk,:)=allD;
FP2(1:nj,1:nk,:)=FPrat;
DP2(1:nj,1:nk,:)=DPrat;
Tgge2(1:nj,1:nk,:)=Tgge;
Tflev2(1:nj,1:nk,:)=Tflev;

%%
% colors
cmBP=cbrewer('seq','BuPu',50);

for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    
    %% Sum mean biom over stages
    f2=figure(2);
    subplot(4,4,s)
    pcolor(agrid,jgrid,squeeze(log10(allF2(:,:,s))))
    colorbar
    colormap(cmBP)
    caxis([-15 0])
    %set(gca,'XTick',mets,'XTickLabel',mets,...
    %    'YTick',cmaxs,'YTickLabel',cmaxs)
    if (s==2)
        title({'log10 Mean F Biom (g m^-^2) in final year'; loc})
    else
        title(loc)
    end
    xlabel('acmax')
    ylabel('aenc')
    stamp(pname)
    
    f3=figure(3);
    subplot(4,4,s)
    pcolor(agrid,jgrid,squeeze(log10(allP2(:,:,s))))
    colorbar
    colormap(cmBP)
    caxis([-2 2])
    %     set(gca,'XTick',mets,'XTickLabel',mets,...
    %         'YTick',cmaxs,'YTickLabel',cmaxs)
    if (s==2)
        title({'log10 Mean P Biom (g m^-^2) in final year'; loc})
    else
        title(loc)
    end
    xlabel('acmax')
    ylabel('aenc')
    stamp(pname)
    
    
    f4=figure(4);
    subplot(4,4,s)
    pcolor(agrid,jgrid,squeeze(log10(allD2(:,:,s))))
    colorbar
    colormap(cmBP)
    caxis([-2 2])
    %     set(gca,'XTick',mets,'XTickLabel',mets,...
    %         'YTick',cmaxs,'YTickLabel',cmaxs)
    if (s==2)
        title({'log10 Mean D Biom (g m^-^2) in final year'; loc})
    else
        title(loc)
    end
    xlabel('acmax')
    ylabel('aenc')
    stamp(pname)
    
    f5=figure(5);
    subplot(4,4,s)
    pcolor(agrid,jgrid,squeeze(FP2(:,:,s)))
    colorbar
    cmocean('balance')
    caxis([0 1])
    %     set(gca,'XTick',mets,'XTickLabel',mets,...
    %         'YTick',cmaxs,'YTickLabel',cmaxs)
    if (s==2)
        title({'Frac F:P in final year'; loc})
    else
        title(loc)
    end
    xlabel('acmax')
    ylabel('aenc')
    stamp(pname)
    
    f6=figure(6);
    subplot(4,4,s)
    pcolor(agrid,jgrid,squeeze(DP2(:,:,s)))
    colorbar
    cmocean('balance')
    caxis([0 1])
    %     set(gca,'XTick',mets,'XTickLabel',mets,...
    %         'YTick',cmaxs,'YTickLabel',cmaxs)
    if (s==2)
        title({'Frac D:P in final year'; loc})
    else
        title(loc)
    end
    xlabel('acmax')
    ylabel('aenc')
    stamp(pname)
    
    %% Feeding level
    f7=figure(7);
    subplot(4,4,s)
    pcolor(agrid,jgrid,squeeze(Tflev2(:,:,s)))
    colorbar
    colormap(cmBP)
    caxis([0.5 1])
    %set(gca,'XTick',mets,'XTickLabel',mets,...
    %    'YTick',cmaxs,'YTickLabel',cmaxs)
    if (s==2)
        title({'Mean feeding level'; loc})
    else
        title(loc)
    end
    xlabel('acmax')
    ylabel('aenc')
    stamp(pname)
    
    %% GGE
    f8=figure(8);
    subplot(4,4,s)
    pcolor(agrid,jgrid,squeeze(Tgge2(:,:,s)))
    colorbar
    colormap(cmBP)
    caxis([0 1])
    %set(gca,'XTick',mets,'XTickLabel',mets,...
    %    'YTick',cmaxs,'YTickLabel',cmaxs)
    if (s==2)
        title({'Mean gross growth efficiency'; loc})
    else
        title(loc)
    end
    xlabel('acmax')
    ylabel('aenc')
    stamp(pname)
    
end %spots
%%
print(f2,'-dpng',[figp pname '_totF_mean_biomass_type_all_locs.png'])
print(f3,'-dpng',[figp pname '_totP_mean_biomass_type_all_locs.png'])
print(f4,'-dpng',[figp pname '_totD_mean_biomass_type_all_locs.png'])
print(f5,'-dpng',[figp pname '_FP_frac_all_locs.png'])
print(f6,'-dpng',[figp pname '_DP_frac_all_locs.png'])
print(f7,'-dpng',[figp pname '_flev_all_locs.png'])
print(f8,'-dpng',[figp pname '_gge_all_locs.png'])



