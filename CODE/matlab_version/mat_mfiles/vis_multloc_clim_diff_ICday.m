% Visualize output of POEM IC day tests
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as mat files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_area_1deg.mat']);
AREA_OCN = max(area,1);

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
load('MyColormaps.mat')
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

%% Loop over kts
% POEM file info
icd = 15:45:365;
Pfish = NaN*ones(length(ID),length(icd)+1);
Ffish = Pfish;
Dfish = Pfish;
for i=1:length(icd)
    dayic = icd(i);
    
    sfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,...
        '/Climatol_All_fish03_ICday_',num2str(dayic),'.mat'];
    load(sfile,...
        'sf_mean','mf_mean','sp_mean','mp_mean','lp_mean',...
        'sd_mean','md_mean','ld_mean');
    
    Ffish(:,i) = sf_mean+mf_mean;
    Pfish(:,i) = sp_mean+mp_mean+lp_mean;
    Dfish(:,i) = sd_mean+md_mean+ld_mean;
    
    clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean
end

%% Base to compare to
sfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,...
    '/Means_bio_prod_fish_Climatol_All_fish03_',cfile,'.mat'];
load(sfile,...
    'sf_mean','mf_mean','sp_mean','mp_mean','lp_mean',...
    'sd_mean','md_mean','ld_mean');

Ffish(:,i+1) = sf_mean+mf_mean;
Pfish(:,i+1) = sp_mean+mp_mean+lp_mean;
Dfish(:,i+1) = sd_mean+md_mean+ld_mean;

clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

%% Start with all fish
Afish = Ffish+Pfish+Dfish;

A1=NaN*ones(ni,nj);
A2=NaN*ones(ni,nj);
A3=NaN*ones(ni,nj);
A4=NaN*ones(ni,nj);
A5=NaN*ones(ni,nj);
A6=NaN*ones(ni,nj);
A7=NaN*ones(ni,nj);
A8=NaN*ones(ni,nj);
A9=NaN*ones(ni,nj);

A1(ID)=Afish(:,1);
A2(ID)=Afish(:,2);
A3(ID)=Afish(:,3);
A4(ID)=Afish(:,4);
A5(ID)=Afish(:,5);
A6(ID)=Afish(:,6);
A7(ID)=Afish(:,7);
A8(ID)=Afish(:,8);
A9(ID)=Afish(:,9);

diff1 = (A1-A9)./A9;
diff2 = (A2-A9)./A9;
diff3 = (A3-A9)./A9;
diff4 = (A4-A9)./A9;
diff5 = (A5-A9)./A9;
diff6 = (A6-A9)./A9;
diff7 = (A7-A9)./A9;
diff8 = (A8-A9)./A9;

%% Maps
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diff1)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.25 0.25]);
colorbar
set(gcf,'renderer','painters')
title(['IC day ' num2str(icd(1))]);
print('-dpng',[pp 'Climatol_' harv '_ICday_comp_',num2str(icd(1)),'.png'])

%2
figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diff2)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.25 0.25]);
colorbar
set(gcf,'renderer','painters')
title(['IC day ' num2str(icd(2))]);
print('-dpng',[pp 'Climatol_' harv '_ICday_comp_',num2str(icd(2)),'.png'])

%3
figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diff3)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.25 0.25]);
colorbar
set(gcf,'renderer','painters')
title(['IC day ' num2str(icd(3))]);
print('-dpng',[pp 'Climatol_' harv '_ICday_comp_',num2str(icd(3)),'.png'])

%4
figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diff4)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.25 0.25]);
colorbar
set(gcf,'renderer','painters')
title(['IC day ' num2str(icd(4))]);
print('-dpng',[pp 'Climatol_' harv '_ICday_comp_',num2str(icd(4)),'.png'])

%5
figure(5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diff5)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.25 0.25]);
colorbar
set(gcf,'renderer','painters')
title(['IC day ' num2str(icd(5))]);
print('-dpng',[pp 'Climatol_' harv '_ICday_comp_',num2str(icd(5)),'.png'])

%6
figure(6)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diff6)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.25 0.25]);
colorbar
set(gcf,'renderer','painters')
title(['IC day ' num2str(icd(6))]);
print('-dpng',[pp 'Climatol_' harv '_ICday_comp_',num2str(icd(6)),'.png'])

figure(7)
%7
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diff7)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.25 0.25]);
colorbar
set(gcf,'renderer','painters')
title(['IC day ' num2str(icd(7))]);
print('-dpng',[pp 'Climatol_' harv '_ICday_comp_',num2str(icd(7)),'.png'])

%8
figure(8)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diff8)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.25 0.25]);
set(gcf,'renderer','painters')
title(['IC day ' num2str(icd(8))]);
colorbar
print('-dpng',[pp 'Climatol_' harv '_ICday_comp_',num2str(icd(8)),'.png'])

%% Calc differences in total biomass
A1tot = A1 .* AREA_OCN;
A2tot = A2 .* AREA_OCN;
A3tot = A3 .* AREA_OCN;
A4tot = A4 .* AREA_OCN;
A5tot = A5 .* AREA_OCN;
A6tot = A6 .* AREA_OCN;
A7tot = A7 .* AREA_OCN;
A8tot = A8 .* AREA_OCN;
A9tot = A9 .* AREA_OCN;

tot_bioA1 = nansum(A1(:));
tot_bioA2 = nansum(A2(:));
tot_bioA3 = nansum(A3(:));
tot_bioA4 = nansum(A4(:));
tot_bioA5 = nansum(A5(:));
tot_bioA6 = nansum(A6(:));
tot_bioA7 = nansum(A7(:));
tot_bioA8 = nansum(A8(:));
tot_bioA9 = nansum(A9(:));

dbiom(:,1) = icd';
dbiom(1,2) = (tot_bioA1 - tot_bioA9) ./ tot_bioA9;
dbiom(2,2) = (tot_bioA2 - tot_bioA9) ./ tot_bioA9;
dbiom(3,2) = (tot_bioA3 - tot_bioA9) ./ tot_bioA9;
dbiom(4,2) = (tot_bioA4 - tot_bioA9) ./ tot_bioA9;
dbiom(5,2) = (tot_bioA5 - tot_bioA9) ./ tot_bioA9;
dbiom(6,2) = (tot_bioA6 - tot_bioA9) ./ tot_bioA9;
dbiom(7,2) = (tot_bioA7 - tot_bioA9) ./ tot_bioA9;
dbiom(8,2) = (tot_bioA8 - tot_bioA9) ./ tot_bioA9;

dbiom(1,3) = nanmax(diff1(:));
dbiom(2,3) = nanmax(diff2(:));
dbiom(3,3) = nanmax(diff3(:));
dbiom(4,3) = nanmax(diff4(:));
dbiom(5,3) = nanmax(diff5(:));
dbiom(6,3) = nanmax(diff6(:));
dbiom(7,3) = nanmax(diff7(:));
dbiom(8,3) = nanmax(diff8(:));

dbiom(1,4) = nanmin(diff1(:));
dbiom(2,4) = nanmin(diff2(:));
dbiom(3,4) = nanmin(diff3(:));
dbiom(4,4) = nanmin(diff4(:));
dbiom(5,4) = nanmin(diff5(:));
dbiom(6,4) = nanmin(diff6(:));
dbiom(7,4) = nanmin(diff7(:));
dbiom(8,4) = nanmin(diff8(:));

dbiom(1,5) = nanmean(diff1(:));
dbiom(2,5) = nanmean(diff2(:));
dbiom(3,5) = nanmean(diff3(:));
dbiom(4,5) = nanmean(diff4(:));
dbiom(5,5) = nanmean(diff5(:));
dbiom(6,5) = nanmean(diff6(:));
dbiom(7,5) = nanmean(diff7(:));
dbiom(8,5) = nanmean(diff8(:));

Tab=array2table(dbiom,'VariableNames',{'ICday','TotalDiffg','MaxDiffgm2',...
    'MinDiffgm2','MeanDiffgm2'});
writetable(Tab,[fpath 'Climatol_' harv '_ICday_comp.csv'],'Delimiter',',','WriteRowNames',true);
save([fpath 'Climatol_' harv '_ICday_comp.mat']);







