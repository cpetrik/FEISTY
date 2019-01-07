% Visualize output of POEM biol rate eq tests
% Spinup at one location
% 150 years, monthly means saved
% Independently change coeffs for cmax and met fns

clear all
close all

%GFDL/NC/Matlab_new_size/Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100/param_sens/
cfile = 'Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100';
   
pp = ['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];
figp = [pp cfile '/param_sens/'];

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
spots=spots';

nfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,'/param_sens/'];
load([nfile 'Locs_Climatol_All_fish03_means_AJsearch.mat']);
pname = 'Climatol_All_fish03_means_AJsearch';


%%
aparam = 0.5:0.1:1;
jparam = 0.5:0.1:1;


% Biomass of each type
allF = SF+MF;
allP = SP+MP+LP;
allD = SD+MD+LD;
allB = BI;

%%
FPrat = squeeze(allF./(allF+allP));
DPrat = squeeze(allD./(allD+allP));

jays = [jparam 1.1];
ays = [aparam 1.1];
[jgrid,agrid]=meshgrid(jays,ays);

allF2 = NaN*ones(7,7,16);
allP2 = NaN*ones(7,7,16);
allD2 = NaN*ones(7,7,16);
FP2 = NaN*ones(7,7,16);
DP2 = NaN*ones(7,7,16);

allF2(1:6,1:6,:)=allF;
allP2(1:6,1:6,:)=allP;
allD2(1:6,1:6,:)=allD;
FP2(1:6,1:6,:)=FPrat;
DP2(1:6,1:6,:)=DPrat;

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
    xlabel('A pref')
    ylabel('J pref')
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
    xlabel('A pref')
    ylabel('J pref')
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
    xlabel('A pref')
    ylabel('J pref')
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
    xlabel('A pref')
    ylabel('J pref')
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
    xlabel('A pref')
    ylabel('J pref')
    stamp(pname)
    
end %spots
%%
print(f2,'-dpng',[figp pname '_totF_mean_biomass_type_all_locs.png'])
print(f3,'-dpng',[figp pname '_totP_mean_biomass_type_all_locs.png'])
print(f4,'-dpng',[figp pname '_totD_mean_biomass_type_all_locs.png'])
print(f5,'-dpng',[figp pname '_FP_frac_all_locs.png'])
print(f6,'-dpng',[figp pname '_DP_frac_all_locs.png'])

