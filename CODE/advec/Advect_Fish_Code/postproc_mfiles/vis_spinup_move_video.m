% Visualize output of Spinup Y1
% Saved as mat files

clear 
close all

%%
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/CORE/'];
ppath = [pp cfile '/CORE/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%%
mod = 'Spinup1988_move_prey_v21_dt24h_All_fish03_Y1';
load([fpath mod '.mat']);
%load([fpath 'Means_' exper cfile '.mat']);

save([fpath mod '.mat'],'S_Bent_bio','S_Sml_f','S_Sml_p','S_Sml_d','S_Med_f',...
    'S_Med_p','S_Med_d','S_Lrg_p','S_Lrg_d','GRD1','GRD2','exper');

%%
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

[ni,nj]=size(geolon_t);
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; 
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% colors
cm10=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black
    

set(groot,'defaultAxesColorOrder',cm10);

%% Take means
time = 1:365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];


%% video
%whole yr
Vmf=NaN*ones(ni,nj,37);
vmf=NaN*ones(ni,nj);
yr = 1:10:365;
for t= 1:length(yr)
    yd = yr(t);
    vmf=NaN*ones(ni,nj);
    vmf(GRD.ID) = (S_Med_f(:,yd));
    Vmf(:,:,t) = vmf;
end

%%
figure(20)
h = imagescn(geolat_t',geolon_t',log10(Vmf(:,:,1))');
cb = colorbar;
ylabel(cb,'log10 abund (g m^-^2)')
cmocean dense
title('Forage fish')
%title(datestr(datenum(0,0,1),'dd'))
clim([-2 2])

hold on
% he = earthimage;
% uistack(he,'bottom')

gif([ppath mod 'MF_Y1.gif'],'frame',gcf,'delaytime',1/12,'nodither')

for k=2:37
    yd = yr(k);
    h.CData = real(log10(Vmf(:,:,k))');
    %title(datestr(datenum(0,0,yd),'dd'))
    gif
end

%%
Vlp=NaN*ones(ni,nj);
Vld=NaN*ones(ni,nj);

Vlp(GRD.ID)=mean(S_Lrg_p(:,1),2,'omitnan');
Vld(GRD.ID)=mean(S_Lrg_d(:,1),2,'omitnan');
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Vmf))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 mean F (g m^-^2) Day 1')

for t=2:365
    Vmf=NaN*ones(ni,nj);
    Vlp=NaN*ones(ni,nj);
    Vld=NaN*ones(ni,nj);

    Vmf(GRD.ID)=mean(S_Med_f(:,t),2,'omitnan');
    Vlp(GRD.ID)=mean(S_Lrg_p(:,t),2,'omitnan');
    Vld(GRD.ID)=mean(S_Lrg_d(:,t),2,'omitnan');



end


