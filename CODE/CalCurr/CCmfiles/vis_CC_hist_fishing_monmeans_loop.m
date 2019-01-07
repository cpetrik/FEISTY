% Visualize output of POEM
% ESM2M Historic (1861-2005) at all locations
% 145 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

kays = 0.0405:0.01:0.0916;
bees = 0.125:0.025:0.25; %bees = 0.1:0.05:0.35;
harv = '03';

for g = 1:length(bees)
    bpow = bees(g);
    
    for n = 1:length(kays)
        kt = kays(n);
        
        tkfn = num2str(100+int64(100*kt));
        tbfn = num2str(100+int64(100*bpow));
        
        cfile = ['Dc_enc70_cmax-metab20_b',tbfn(2:end),'_k',tkfn(2:end),'_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100'];
        
        fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
        ppath = [pp cfile '/'];
        if (~isdir(ppath))
            mkdir(ppath)
        end
        
        load([fpath 'Means_Hindcast_fished_' harv '_' cfile '.mat']);
        
        close all
        
        %% colors
        cm9=[0.5 0.5 0;... %tan/army
            0 0.7 0;...   %g
            1 0 1;...     %m
            1 0 0;...     %r
            0.5 0 0;...   %maroon
            0/255 206/255 209/255;... %turq
            0 0.5 0.75;...   %med blue
            0 0 0.75;...  %b
            0 0 0];...      %black
            
        cm21=[1 0.5 0;...   %orange
            0.5 0.5 0;... %tan/army
            0 0.7 0;...   %g
            0 1 1;...     %c
            0 0 0.75;...  %b
            0.5 0 1;...   %purple
            1 0 1;...     %m
            1 0 0;...     %r
            0.5 0 0;...   %maroon
            0.75 0.75 0.75;... %lt grey
            0.5 0.5 0.5;...    %med grey
            49/255 79/255 79/255;... %dk grey
            0 0 0;...      %black
            1 1 0;...      %yellow
            127/255 255/255 0;... %lime green
            0 0.5 0;...    %dk green
            0/255 206/255 209/255;... %turq
            0 0.5 0.75;...   %med blue
            188/255 143/255 143/255;... %rosy brown
            255/255 192/255 203/255;... %pink
            255/255 160/255 122/255]; %peach
        
        set(groot,'defaultAxesColorOrder',cm9);
        
        %% Plots in time
        y = time;
        nt = length(time);
        
        % % Piscivore
        % figure(1)
        % subplot(4,1,1)
        % plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
        % plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
        % plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
        % xlim([y(1) y(end)])
        % title(['Pelagic Piscivores Climatol ' tharv])
        % ylabel('log10 Biomass (g m^-^2)')
        % legend('Larvae','Juveniles','Adults')
        % legend('location','southeast')
        % %stamp(cfile)
        %
        % subplot(4,1,2)
        % plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
        % xlim([y(1) y(end)])
        % title('Larvae')
        % ylabel('log10 Biomass (g m^-^2)')
        %
        % subplot(4,1,3)
        % plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
        % xlim([y(1) y(end)])
        % title('Juveniles')
        % ylabel('log10 Biomass (g m^-^2)')
        %
        % subplot(4,1,4)
        % plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
        % xlim([y(1) y(end)])
        % title('Adults')
        % xlabel('Time (mo)')
        % ylabel('log10 Biomass (g m^-^2)')
        % print('-dpng',[ppath 'Hindcast_fished_' harv '_P_time.png'])
        %
        % % Planktivore
        % sf_tmean=sf_tmean(1:length(y));
        % figure(2)
        % subplot(3,1,1)
        % plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
        % plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
        % xlim([y(1) y(end)])
        % title(['Forage Fishes Climatol ' tharv])
        % xlabel('Time (mo)')
        % ylabel('log10 Biomass (g m^-^2)')
        % legend('Immature','Adults')
        % legend('location','southeast')
        % %stamp(cfile)
        %
        % subplot(3,1,2)
        % plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
        % xlim([y(1) y(end)])
        % title('Immature')
        % ylabel('log10 Biomass (g m^-^2)')
        %
        % subplot(3,1,3)
        % plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
        % xlim([y(1) y(end)])
        % title('Adults')
        % xlabel('Time (mo)')
        % ylabel('log10 Biomass (g m^-^2)')
        % print('-dpng',[ppath 'Hindcast_fished_' harv '_F_time.png'])
        %
        % % Detritivore
        % figure(3)
        % subplot(4,1,1)
        % plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
        % plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
        % plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
        % xlim([y(1) y(end)])
        % title(['Demersal Piscivores Climatol ' tharv])
        % ylabel('log10 Biomass (g m^-^2)')
        % legend('Larvae','Juveniles','Adults')
        % legend('location','southeast')
        % %stamp(cfile)
        %
        % subplot(4,1,2)
        % plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
        % xlim([y(1) y(end)])
        % title('Larvae')
        % ylabel('log10 Biomass (g m^-^2)')
        %
        % subplot(4,1,3)
        % plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
        % xlim([y(1) y(end)])
        % title('Juveniles')
        % ylabel('log10 Biomass (g m^-^2)')
        %
        % subplot(4,1,4)
        % plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
        % xlim([y(1) y(end)])
        % title('Adults')
        % xlabel('Time (mo)')
        % ylabel('log10 Biomass (g m^-^2)')
        % print('-dpng',[ppath 'Hindcast_fished_' harv '_D_time.png'])
        
        %% All size classes of all
        
        figure(4)
        plot(y,log10(sf_tmean(1:nt)),'Linewidth',1); hold on;
        plot(y,log10(mf_tmean),'Linewidth',1); hold on;
        plot(y,log10(sp_tmean),'Linewidth',1); hold on;
        plot(y,log10(mp_tmean),'Linewidth',1); hold on;
        plot(y,log10(lp_tmean),'Linewidth',1); hold on;
        plot(y,log10(sd_tmean),'Linewidth',1); hold on;
        plot(y,log10(md_tmean),'Linewidth',1); hold on;
        plot(y,log10(ld_tmean),'Linewidth',1); hold on;
        legend('SF','MF','SP','MP','LP','SD','MD','LD')
        legend('location','eastoutside')
        xlim([y(1) y(end)])
        ylim([-5 2])
        xlabel('Time (mo)')
        ylabel('log10 Biomass (g m^-^2)')
        title(['Hindcast fishing ' harv])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished_' harv '_all_sizes.png'])
        
        figure(5)
        F = sf_tmean(1:nt)+mf_tmean;
        P = sp_tmean+mp_tmean+lp_tmean;
        D = sd_tmean+md_tmean+ld_tmean;
        
        plot(y,log10(F),'r','Linewidth',2); hold on;
        plot(y,log10(P),'b','Linewidth',2); hold on;
        plot(y,log10(D),'k','Linewidth',2); hold on;
        legend('F','P','D')
        legend('location','eastoutside')
        xlim([y(1) y(end)])
        ylim([-5 2])
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        title(['Hindcast fishing ' harv])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished_' harv '_all_types.png'])
        
        %% FISHING All size classes of all
        
        figure(6)
        plot(y,log10(mf_tmy),'color',[0 0.7 0],'Linewidth',1); hold on;
        plot(y,log10(mp_tmy),'color',[1 0 0],'Linewidth',1); hold on;
        plot(y,log10(lp_tmy),'color',[0.5 0 0],'Linewidth',1); hold on;
        plot(y,log10(md_tmy),'color',[0 0.5 0.75],'Linewidth',1); hold on;
        plot(y,log10(ld_tmy),'color',[0 0 0.75],'Linewidth',1); hold on;
        legend('MF','MP','LP','MD','LD')
        legend('location','eastoutside')
        xlim([y(1) y(end)])
        ylim([-7 0])
        xlabel('Time (mo)')
        ylabel('log10 Catch (g m^-^2)')
        title(['Hindcast fishing ' harv])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished_' harv '_catch_all_sizes.png'])
        
        figure(7)
        F = mf_tmy;
        P = mp_tmy+lp_tmy;
        D = md_tmy+ld_tmy;
        
        plot(y,log10(F),'r','Linewidth',2); hold on;
        plot(y,log10(P),'b','Linewidth',2); hold on;
        plot(y,log10(D),'k','Linewidth',2); hold on;
        legend('F','P','D')
        legend('location','eastoutside')
        xlim([y(1) y(end)])
        ylim([-7 0])
        xlabel('Time (y)')
        ylabel('log10 Catch (g m^-^2)')
        title(['Hindcast fishing ' harv])
        print('-dpng',[ppath 'Hindcast_fished_' harv '_catch_all_types.png'])
        
        
        %% Plots in space
        [ni,nj]=size(geolon_t);
        
        Zsf=NaN*ones(ni,nj);
        Zsp=NaN*ones(ni,nj);
        Zsd=NaN*ones(ni,nj);
        Zmf=NaN*ones(ni,nj);
        Zmp=NaN*ones(ni,nj);
        Zmd=NaN*ones(ni,nj);
        Zlp=NaN*ones(ni,nj);
        Zld=NaN*ones(ni,nj);
        Zb=NaN*ones(ni,nj);
        
        Zsf(grid(:,1))=sf_mean;
        Zsp(grid(:,1))=sp_mean;
        Zsd(grid(:,1))=sd_mean;
        Zmf(grid(:,1))=mf_mean;
        Zmp(grid(:,1))=mp_mean;
        Zmd(grid(:,1))=md_mean;
        Zlp(grid(:,1))=lp_mean;
        Zld(grid(:,1))=ld_mean;
        Zb(grid(:,1))=b_mean;
        
        Cmf=NaN*ones(ni,nj);
        Cmp=NaN*ones(ni,nj);
        Cmd=NaN*ones(ni,nj);
        Clp=NaN*ones(ni,nj);
        Cld=NaN*ones(ni,nj);
        
        Cmf(grid(:,1))=mf_my;
        Cmp(grid(:,1))=mp_my;
        Cmd(grid(:,1))=md_my;
        Clp(grid(:,1))=lp_my;
        Cld(grid(:,1))=ld_my;
        
        %% bent
        figure(9)
        surf(geolon_t,geolat_t,log10(Zb)); view(2); hold on;
        shading flat
        title('Historic fished 1956-2005 log10 mean benthic biomass (g m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([-2.5 0.5])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished',harv,'_global_BENT.png'])
        
        %
        mgZb = (Zb/9)*1e3;
        figure(10)
        surf(geolon_t,geolat_t,log10(mgZb)); view(2); hold on;
        shading flat
        title('Historic fished 1956-2005 log10 mean benthic biomass (mg C m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([-0.8 2.3])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished',harv,'_global_BENT_mgC.png'])
        
%         % sp
%         figure(11)
%         surf(geolon_t,geolat_t,log10(Zsp)); view(2); hold on;
%         shading flat
%         title('Historic fished 1956-2005 log10 mean Larval P biomass (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-2 1])
%         stamp(cfile)
%         print('-dpng',[ppath 'Hindcast_fished',harv,'_global_SP.png'])
%         
%         % sf
%         figure(12)
%         surf(geolon_t,geolat_t,log10(Zsf)); view(2); hold on;
%         shading flat
%         title('Historic fished 1956-2005 log10 mean Larval F biomass (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-2 1])
%         stamp(cfile)
%         print('-dpng',[ppath 'Hindcast_fished',harv,'_global_SF.png'])
%         
%         % sd
%         figure(13)
%         surf(geolon_t,geolat_t,log10(Zsd)); view(2); hold on;
%         shading flat
%         title('Historic fished 1956-2005 log10 mean Larval D biomass (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-2 1])
%         stamp(cfile)
%         print('-dpng',[ppath 'Hindcast_fished',harv,'_global_SD.png'])
%         
%         % mp
%         figure(14)
%         surf(geolon_t,geolat_t,log10(Zmp)); view(2); hold on;
%         shading flat
%         title('Historic fished 1956-2005 log10 mean Juvenile P biomass (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-2 1])
%         stamp(cfile)
%         print('-dpng',[ppath 'Hindcast_fished',harv,'_global_MP.png'])
%         
%         % mf
%         figure(15)
%         surf(geolon_t,geolat_t,log10(Zmf)); view(2); hold on;
%         shading flat
%         title('Historic fished 1956-2005 log10 mean Adult F biomass (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-2 1])
%         stamp(cfile)
%         print('-dpng',[ppath 'Hindcast_fished',harv,'_global_MF.png'])
%         
%         % md
%         figure(16)
%         surf(geolon_t,geolat_t,log10(Zmd)); view(2); hold on;
%         shading flat
%         title('Historic fished 1956-2005 log10 mean Juvenile D biomass (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-2 1])
%         stamp(cfile)
%         print('-dpng',[ppath 'Hindcast_fished',harv,'_global_MD.png'])
%         
%         % lp
%         figure(17)
%         surf(geolon_t,geolat_t,log10(Zlp)); view(2); hold on;
%         shading flat
%         title('Historic fished 1956-2005 log10 mean Adult P biomass (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-2 1])
%         stamp(cfile)
%         print('-dpng',[ppath 'Hindcast_fished',harv,'_global_LP.png'])
%         
%         % ld
%         figure(18)
%         surf(geolon_t,geolat_t,log10(Zld)); view(2); hold on;
%         shading flat
%         title('Historic fished 1956-2005 log10 mean Adult D biomass (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-2 1])
%         stamp(cfile)
%         print('-dpng',[ppath 'Hindcast_fished',harv,'_global_LD.png'])
        
        %% Diff maps of all fish
        All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
        AllF = Zsf+Zmf;
        AllP = Zsp+Zmp+Zlp;
        AllD = Zsd+Zmd+Zld;
        AllS = Zsp+Zsf+Zsd;
        AllM = Zmp+Zmf+Zmd;
        AllL = Zlp+Zld;
        FracPD = AllP ./ (AllP+AllD);
        FracPF = AllP ./ (AllP+AllF);
        FracPFvD = (AllP+AllF) ./ (AllP+AllF+AllD);
        
        %% ALL
        figure(21)
        surf(geolon_t,geolat_t,log10(All)); view(2); hold on;
        shading flat
        title('Historic fished 1956-2005 log10 mean biomass All Fishes (g m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([-1 1])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished',harv,'_global_All.png'])
        
        % all F
        figure(22)
        surf(geolon_t,geolat_t,log10(AllF)); view(2); hold on;
        shading flat
        title('Historic fished 1956-2005 log10 mean biomass All F (g m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([-2 1])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished',harv,'_global_AllF.png'])
        
        % all D
        figure(23)
        surf(geolon_t,geolat_t,log10(AllD)); view(2); hold on;
        shading flat
        title('Historic fished 1956-2005 log10 mean biomass All D (g m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([-2 1])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished',harv,'_global_AllD.png'])
        
        % All P
        figure(24)
        surf(geolon_t,geolat_t,log10(AllP)); view(2); hold on;
        shading flat
        title('Historic fished 1956-2005 log10 mean biomass All P (g m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([-2 1])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished',harv,'_global_AllP.png'])
        
        % FracPD
        figure(25)
        surf(geolon_t,geolat_t,FracPD); view(2); hold on;
        shading flat
        title('Historic fished 1956-2005 P:D mean biomass(g m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([0 1])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished',harv,'_global_FracPD.png'])
        
        % FracPF
        figure(26)
        surf(geolon_t,geolat_t,FracPF); view(2); hold on;
        shading flat
        title('Historic fished 1956-2005 P:F mean biomass (g m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([0 1])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished',harv,'_global_FracPF.png'])
        
        
        %% CATCH
%         % mp
%         figure(43)
%         surf(geolon_t,geolat_t,log10(Cmp)); view(2); hold on;
%         shading flat
%         title('log10 mean Juvenile P catch (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-5 -2])
%         %stamp(cfile)
%         print('-dpng',[ppath 'Hindcast_fished',harv, '_global_MP_catch.png'])
%         
%         % mf
%         figure(44)
%         surf(geolon_t,geolat_t,log10(Cmf)); view(2); hold on;
%         shading flat
%         title('log10 mean Adult F catch (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-5 -2])
%         %stamp(cfile)
%         print('-dpng',[ppath 'Hindcast_fished',harv, '_global_MF_catch.png'])
%         
%         % md
%         figure(45)
%         surf(geolon_t,geolat_t,log10(Cmd)); view(2); hold on;
%         shading flat
%         title('log10 mean Juvenile D catch (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-5 -2])
%         %stamp(cfile)
%         print('-dpng',[ppath 'Hindcast_fished',harv, '_global_MD_catch.png'])
%         
%         % lp
%         figure(46)
%         surf(geolon_t,geolat_t,log10(Clp)); view(2); hold on;
%         shading flat
%         title('log10 mean Adult P catch (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-5 -2])
%         %stamp(cfile)
%         print('-dpng',[ppath 'Hindcast_fished',harv, '_global_LP_catch.png'])
%         
%         % ld
%         figure(47)
%         surf(geolon_t,geolat_t,log10(Cld)); view(2); hold on;
%         shading flat
%         title('log10 mean Adult D catch (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-5 -2])
%         %stamp(cfile)
%         print('-dpng',[ppath 'Hindcast_fished',harv, '_global_LD_catch.png'])
%         
        %% Diff maps of all fish
        CAll = Cmp+Cmf+Cmd+Clp+Cld;
        CAllF = Cmf;
        CAllP = Cmp+Clp;
        CAllD = Cmd+Cld;
        CAllM = Cmp+Cmf+Cmd;
        CAllL = Clp+Cld;
        
        % ALL
        figure(48)
        surf(geolon_t,geolat_t,log10(CAll)); view(2); hold on;
        shading flat
        title('log10 mean catch All Fishes (g m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([-5 -2])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished',harv, '_global_All_catch.png'])
        
        % all F
        figure(49)
        surf(geolon_t,geolat_t,log10(CAllF)); view(2); hold on;
        shading flat
        title('log10 mean catch All F (g m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([-5 -2])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished',harv, '_global_AllF_catch.png'])
        
        % all D
        figure(50)
        surf(geolon_t,geolat_t,log10(CAllD)); view(2); hold on;
        shading flat
        title('log10 mean catch All D (g m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([-5 -2])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished',harv, '_global_AllD_catch.png'])
        
        % All P
        figure(51)
        surf(geolon_t,geolat_t,log10(CAllP)); view(2); hold on;
        shading flat
        title('log10 mean catch All P (g m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([-5 -2])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished',harv, '_global_AllP_catch.png'])
        
        % all M
        figure(52)
        surf(geolon_t,geolat_t,log10(CAllM)); view(2); hold on;
        shading flat
        title('log10 mean catch All M (g m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([-5 -2])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished',harv, '_global_AllM_catch.png'])
        
        % All L
        figure(53)
        surf(geolon_t,geolat_t,log10(CAllL)); view(2); hold on;
        shading flat
        title('log10 mean catch All L (g m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([-5 -2])
        stamp(cfile)
        print('-dpng',[ppath 'Hindcast_fished',harv, '_global_AllL_catch.png'])
        
    end
end
