% Visualize output of POEM biol rate eq tests
% Spinup at one location
% 150 years, monthly means saved
% Independently change coeffs for cmax and met fns

clear all
close all

%?/GFDL/NC/Matlab_new_size/Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/param_sens/
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
   
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
load([nfile 'Locs_Climatol_All_fish03_means_bm_kt_search.mat']);
pname = 'Climatol_All_fish03_means_bm_kt_search';


%%
bmp = 0.1:0.025:0.325;
ktemp = 0.0605:0.01:0.1205;
% ktemp = 0.0805:0.005:0.1005;
% bees = 0.025:0.025:0.15;
nj=length(bmp);
nk=length(ktemp);

% Biomass of each type
allF = SF+MF;
allP = SP+MP+LP;
allD = SD+MD+LD;
allB = BI;

%%
FPrat = squeeze(allF./(allF+allP));
DPrat = squeeze(allD./(allD+allP));

jays = [bmp 0.35];
ays = [ktemp 0.1305];
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

for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    
    %% Sum mean biom over stages
    f2=figure(2);
    subplot(4,4,s)
    pcolor(agrid,jgrid,squeeze(log10(allF2(:,:,s))))
    colorbar
    colormap(cmBP)
    caxis([-2 2])
    %set(gca,'XTick',mets,'XTickLabel',mets,...
    %    'YTick',cmaxs,'YTickLabel',cmaxs)
    if (s==2)
        title({'log10 Mean F Biom (g m^-^2) in final year'; loc})
    else
        title(loc)
    end
    ylabel('bcmax')
    xlabel('kt')
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
    ylabel('bcmax')
    xlabel('kt')
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
    ylabel('bcmax')
    xlabel('kt')
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
    ylabel('bcmax')
    xlabel('kt')
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
    ylabel('bcmax')
    xlabel('kt')
    stamp(pname)
    
end %spots
%%
print(f2,'-dpng',[figp pname '_totF_mean_biomass_type_all_locs.png'])
print(f3,'-dpng',[figp pname '_totP_mean_biomass_type_all_locs.png'])
print(f4,'-dpng',[figp pname '_totD_mean_biomass_type_all_locs.png'])
print(f5,'-dpng',[figp pname '_FP_frac_all_locs.png'])
print(f6,'-dpng',[figp pname '_DP_frac_all_locs.png'])

