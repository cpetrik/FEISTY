% Visualize output of POEM
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/Bio_rates/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

close all

% plot info
[ni,nj]=size(lon);
geolat_t=lat;
geolon_t=lon;
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));


% colors
cmBP=cbrewer('seq','BuPu',50);

%% Plots in space

cfileA = 'Dc_enc70-b200_m4-b175-k041_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpathA=['/Volumes/GFDL/NC/Matlab_new_size/' cfileA '/'];
load([fpathA 'Climatol_All_fish03.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Afish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Aall=NaN*ones(ni,nj);
Aall(ID)=Afish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

cfileB = 'Dc_enc70-b200_m4-b175-k071_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpathB=['/Volumes/GFDL/NC/Matlab_new_size/' cfileB '/'];
load([fpathB 'Climatol_All_fish03.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Bfish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Ball=NaN*ones(ni,nj);
Ball(ID)=Bfish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

cfileC = 'Dc_enc70-b200_m4-b175-k101_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpathC=['/Volumes/GFDL/NC/Matlab_new_size/' cfileC '/'];
load([fpathC 'Climatol_All_fish03.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Cfish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Call=NaN*ones(ni,nj);
Call(ID)=Cfish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

cfileD = 'Dc_enc70-b200_m4-b175-k131_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpathD=['/Volumes/GFDL/NC/Matlab_new_size/' cfileD '/'];
load([fpathD 'Climatol_All_fish03.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Dfish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Dall=NaN*ones(ni,nj);
Dall(ID)=Dfish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

%% Loop over kts
% POEM file info
harv = 'All_fish03';
kays = [0.0405,0.0555,0.0705,0.0855,0.1005,0.1155,0.1305];
skays={'.04','.0555','.07','.0855','.10','.1155','.13'};

fish = NaN*ones(length(ID),length(kays));

for i=1:length(kays)
    kt=kays(i);
    tkfn = num2str(1000+int64(1000*kt));
    cfile = ['Dc_enc70-b200_m4-b175-k',tkfn(2:end),...
        '_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100'];
    fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
    if (i==4)
        load([fpath 'Means_bio_prod_fish_Climatol_All_fish03_',cfile,'.mat'],...
            'sf_mean','mf_mean','sp_mean',...
            'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
    else
        load([fpath 'Climatol_All_fish03.mat'],...
            'sf_mean','mf_mean','sp_mean',...
            'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
    end

    fish(:,i) = sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+...
        lp_mean+ld_mean;

    clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean
end
%%
A1=NaN*ones(ni,nj);
A2=NaN*ones(ni,nj);
A3=NaN*ones(ni,nj);
A4=NaN*ones(ni,nj);
A5=NaN*ones(ni,nj);
A6=NaN*ones(ni,nj);
A7=NaN*ones(ni,nj);
% A8=NaN*ones(ni,nj);
% A9=NaN*ones(ni,nj);
A1(ID)=fish(:,1);
A2(ID)=fish(:,2);
A3(ID)=fish(:,3);
A4(ID)=fish(:,4);
A5(ID)=fish(:,5);
A6(ID)=fish(:,6);
A7(ID)=fish(:,7);
% A8(ID)=fish(:,8);
% A9(ID)=fish(:,9);

%% 
figure(1)
%A
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(A1))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'A')
text(-0.75,1.75,'kt=0.04')
%B
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(A4))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'B')
text(-0.75,1.75,'kt=0.0855')
%C
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(A7))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'C')
text(-0.75,1.75,'kt=0.13')
%stamp('')
colorbar('Position',[0.04 0.04 0.339 0.025],'orientation','horizontal')
print('-dpng',[pp 'Fig3_Climatol_' harv '_All_comp_kt_4-13_BP.png'])

