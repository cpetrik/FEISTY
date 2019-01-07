% Visualize output of POEM
% Spinup at 100 locations
% 50 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
% pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_Big_sizes/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_OG_sizes/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);


RE = {'1000','0500','0100'};%,'0050','0010','00050','00010'};
reff = [1.0,0.5,0.1];%,0.05,0.01,0.005,0.001];
CarCap = {'050','100','150','200','250','300'};
car = 0.5:0.5:3.0;
benteff = {'05','10'};%,'15','20'};
beff = 0.05:0.05:0.1;
fcrit = 40;
nmort = '0';
kad = 100;
pref = 'D100';

for c = 1:length(CarCap)
    CC = CarCap{c};
    close all
    n=0;
    for R = 1:length(RE)
        rfrac = RE{R};
        for i=1:length(benteff)
            BE = benteff{i};
            
            n=n+1;
            %%
            cfile = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort'  nmort '_BE' BE '_CC' CC '_RE' rfrac];
            cfile2 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort'  nmort '_CC' CC ];
            
            % fpath=['/Volumes/GFDL/NC/Matlab_big_size/' cfile '/'];
            fpath=['/Volumes/GFDL/NC/Matlab_og_size/' cfile '/'];
            ppath = [pp cfile '/'];
            
            load([fpath 'Means_spinup_' cfile '.mat']);
            
            %% Plots in space
            
            [ni,nj]=size(geolon_t);
            
            Zsf=NaN*ones(ni,nj);
            Zsp=NaN*ones(ni,nj);
            Zmf=NaN*ones(ni,nj);
            Zmp=NaN*ones(ni,nj);
            Zlp=NaN*ones(ni,nj);
            
            Zsf(grid(:,1))=sf_mean;
            Zsp(grid(:,1))=sp_mean;
            Zmf(grid(:,1))=mf_mean;
            Zmp(grid(:,1))=mp_mean;
            Zlp(grid(:,1))=lp_mean;
            
            % Diff maps of all fish
            AllF = Zsf+Zmf;
            AllP = Zsp+Zmp+Zlp;
            FracPF = AllP ./ (AllP+AllF);
            
            %% All F
            figure
            m_proj('miller','lat',82);
            surf(geolon_t,geolat_t,log10(AllF)); view(2); hold on;
            shading flat
            % m_coast('patch',[.5 .5 .5],'edgecolor','none');
            % m_grid;
            title('log10 mean biomass All F (g m^-^2)')
            colormap('jet')
            colorbar('h')
            caxis([-2 1])
            xlim([-280 90])
            ylim([-80 90])
            stamp(cfile)
            print('-dpng',[ppath 'Spinup_global_AllF.png'])
            
            % All P
            figure
            m_proj('miller','lat',82);
            surf(geolon_t,geolat_t,log10(AllP)); view(2); hold on;
            shading flat
            % m_coast('patch',[.5 .5 .5],'edgecolor','none');
            % m_grid;
            title('log10 mean biomass All P (g m^-^2)')
            colormap('jet')
            colorbar('h')
            caxis([-2 1])
            xlim([-280 90])
            ylim([-80 90])
            stamp(cfile)
            print('-dpng',[ppath 'Spinup_global_AllP.png'])
            
            % FracPF
            figure
            m_proj('miller','lat',82);
            surf(geolon_t,geolat_t,FracPF); view(2); hold on;
            shading flat
            % m_coast('patch',[.5 .5 .5],'edgecolor','none');
            % m_grid;
            title('P:F mean biomass (g m^-^2)')
            colormap('jet')
            colorbar('h')
            caxis([0 1])
            xlim([-280 90])
            ylim([-80 90])
            stamp(cfile)
            print('-dpng',[ppath 'Spinup_global_FracPF.png'])
            
            % FracPF
            f50 = figure(50);
            subplot(3,2,n)
            surf(geolon_t,geolat_t,FracPF); view(2); hold on;
            shading flat
            title(['P:F BE=' num2str(beff(i))])
            ylabel(['RE=' num2str(reff(R))])
            colormap('jet')
%             colorbar('h')
            caxis([0 1])
            xlim([-280 90])
            ylim([-80 90])
            stamp(cfile)
            
            
        end
    end
    print(f50, '-dpng',[pp cfile2 '_Spinup_global_FracPF.png'])
end


