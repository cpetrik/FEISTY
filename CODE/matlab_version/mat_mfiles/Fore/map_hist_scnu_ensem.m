function map_hist_scnu_ensem(simname,sf_nu50,sp_nu50,sd_nu50,...
    mf_nu50,mp_nu50,md_nu50,lp_nu50,ld_nu50,pp)

% Visualize output of FEISTY

close all

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
grid = csvread([cpath 'grid_csv.csv']);

%% colors
warning off
cmBP=cbrewer('seq','BuPu',50);

%% Plot info
[ni,nj]=size(geolon_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% Plots in space
Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);

Zsf(grid(:,1))=sf_nu50;
Zsp(grid(:,1))=sp_nu50;
Zsd(grid(:,1))=sd_nu50;
Zmf(grid(:,1))=mf_nu50;
Zmp(grid(:,1))=mp_nu50;
Zmd(grid(:,1))=md_nu50;
Zlp(grid(:,1))=lp_nu50;
Zld(grid(:,1))=ld_nu50;

All = (Zmf+Zlp+Zld)/3;
AllF = (Zsf+Zmf)/2;
AllP = (Zsp+Zmp+Zlp)/3;
AllD = (Zsd+Zmd+Zld)/3;
AllS = (Zsp+Zsf+Zsd)/3;
AllM = (Zmp+Zmf+Zmd)/3;
AllL = (Zlp+Zld)/2;

%% 6 plot 
figure(1)
%A 
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Zmf)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.6 0.6]);
colorbar('Position',[0.385 0.695 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Adult Forage','HorizontalAlignment','center')
text(-2.75,1.75,'A')
%B 
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Zlp)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.6 0.6]);
colorbar('Position',[0.385 0.385 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Adult Large pelagic','HorizontalAlignment','center')
text(-2.75,1.75,'B')
% C
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Zld)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.6 0.6]);
colorbar('Position',[0.385 0.075 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Adult Demersal','HorizontalAlignment','center')
text(-2.75,1.75,'C')
% D
subplot('Position',[0.45 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,AllS)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.6 0.6]);
colorbar('Position',[0.825 0.695 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Small','HorizontalAlignment','center')
text(-2.75,1.75,'D')
% E
subplot('Position',[0.45 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,AllM)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.6 0.6]);
colorbar('Position',[0.825 0.385 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Medium','HorizontalAlignment','center')
text(-2.75,1.75,'E')
%F 
subplot('Position',[0.45 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,AllL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.6 0.6]);
colorbar('Position',[0.825 0.075 0.025 0.275],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Large','HorizontalAlignment','center')
text(-2.75,1.75,'F')
print('-dpng',[pp 'Hist_scnu_types_6plot_',simname,'.png'])



end