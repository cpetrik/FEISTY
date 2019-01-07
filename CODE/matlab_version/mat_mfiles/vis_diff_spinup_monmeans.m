Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';

cpath = [Fdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Fdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];
ppath = [pp cfile '/'];

load([Fdrpbx 'Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat'],...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

%% Plots in space
[ni,nj]=size(geolon_t);

st=1:60:length(time);
en=60:60:length(time);
yr=5:5:(length(time)/(12));

for n=1:length(st)
    mf_mean=nanmean(MF.bio(:,st(n):en(n)),2);
    
    Zmf=NaN*ones(ni,nj);
    Zmf(grid(:,1))=mf_mean;
    
    % mf
    figure(15)
    clf
    surf(geolon_t,geolat_t,log10(Zmf)); view(2); hold on;
    shading flat
    title(['Year ' num2str(yr(n)) ' log10 mean Adult F biomass (g m^-^2)'])
    colormap('jet')
    colorbar('h')
    caxis([-2 1])
    print('-dpng',[ppath 'Spinup_global_MF_yr' num2str(yr(n)) '.png'])
end

%% First 10 yrs
st=1:12:120;
en=12:12:120;
yr=1:length(st);

for n=1:length(st)
    mf_mean=nanmean(MF.bio(:,st(n):en(n)),2);
    
    Zmf=NaN*ones(ni,nj);
    Zmf(grid(:,1))=mf_mean;
    
    % mf
    figure(15)
    clf
    surf(geolon_t,geolat_t,log10(Zmf)); view(2); hold on;
    shading flat
    title(['Year ' num2str(yr(n)) ' log10 mean Adult F biomass (g m^-^2)'])
    colormap('jet')
    colorbar('h')
    caxis([-2 1])
    print('-dpng',[ppath 'Spinup_global_MF_yr' num2str(yr(n)) '.png'])
end
