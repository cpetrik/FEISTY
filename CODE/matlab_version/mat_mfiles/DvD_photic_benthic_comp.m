% Z:B ratio by LME 
% Climatology
% 150 years
% Saved as mat files
% Compare to Daniel's model results

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cdir 'temp_100_1deg_ESM26_5yr_clim_191_195.mat'])
load([cdir 'btm_temp_1deg_ESM26_5yr_clim_191_195.mat'])
load([cpath 'LME_clim_temp_zoop_det.mat']);

ptemp_mean_clim=squeeze(nanmean(temp_100,1));
btemp_mean_clim=squeeze(nanmean(btm_temp,1));
tlme = lme_mask_onedeg;
AREA_OCN = max(area,1);

%% POEM
cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';
ppath = [pp cfile '/'];
dpath = [dp cfile '/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat'],'lme_mbio',...
    'lme_mcatch','lme_area','rPD_biom','rPD_catch','rPD_catch_mtkm2');

lme_area_km2 = lme_area * 1e-6;

% POEM LME biomass in MT
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
% MT/km2
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;

plme_rPDcatch = plme_Pmcatch ./ (plme_Pmcatch+plme_Dmcatch);

%% COBALT
% units ESM2.6 in mg C m-2 or mg C m-2 d-1

lme_zl = lme_zl*365;
lme_det = lme_det*365;

FracZDet = lme_z ./ (lme_z+lme_det);
FracZB = lme_z ./ (lme_z+lme_mbio(:,9));
FracZlDet = lme_zl ./ (lme_zl+lme_det);
FracZlB = lme_zl ./ (lme_zl+lme_mbio(:,9));

RatZDet = lme_z ./ (lme_det);
RatZB = lme_z ./ (lme_mbio(:,9));
RatZlDet = lme_zl ./ (lme_det);
RatZlB = lme_zl ./ (lme_mbio(:,9));

%% DvD 
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/DanielVD_PelDem/Colleen_modeledfish_LME.mat')
RatPS = Fphotic ./ Fseabed;

%% Comparison stats
keep=[1:61,63];

%r
rZD=corr(log10(RatPS(keep)),log10(RatZDet(keep)));
rZB=corr(log10(RatPS(keep)),log10(RatZB(keep)));
rZlD=corr(log10(RatPS(keep)),log10(RatZlDet(keep)));
rZlB=corr(log10(RatPS(keep)),log10(RatZlB(keep)));
rZ=corr(log10(Fphotic(keep)),log10(lme_z(keep)));
rZl=corr(log10(Fphotic(keep)),log10(lme_zl(keep)));
rD=corr(log10(Fseabed(keep)),log10(lme_det(keep)));
rB=corr(log10(Fseabed(keep)),log10(lme_mbio(keep,9)));

%root mean square error
o=log10(RatPS(keep));
p=log10(RatZDet(keep));
n = length(o);
num=nansum((p-o).^2);
rmseZD = sqrt(num/n);

o=log10(RatPS(keep));
p=log10(RatZB(keep));
n = length(o);
num=nansum((p-o).^2);
rmseZB = sqrt(num/n);

o=log10(RatPS(keep));
p=log10(RatZlDet(keep));
n = length(o);
num=nansum((p-o).^2);
rmseZlD = sqrt(num/n);

o=log10(RatPS(keep));
p=log10(RatZlB(keep));
n = length(o);
num=nansum((p-o).^2);
rmseZlB = sqrt(num/n);

o=log10(Fphotic(keep));
p=log10(lme_z(keep));
n = length(o);
num=nansum((p-o).^2);
rmseZ = sqrt(num/n);

o=log10(Fphotic(keep));
p=log10(lme_zl(keep));
n = length(o);
num=nansum((p-o).^2);
rmseZl = sqrt(num/n);

o=log10(Fseabed(keep));
p=log10(lme_det(keep));
n = length(o);
num=nansum((p-o).^2);
rmseD = sqrt(num/n);

o=log10(Fseabed(keep));
p=log10(lme_mbio(keep,9));
n = length(o);
num=nansum((p-o).^2);
rmseB = sqrt(num/n);

% Table
fish_stat(1,1) = rZD;
fish_stat(1,2) = rmseZD;
fish_stat(2,1) = rZB;
fish_stat(2,2) = rmseZB;
fish_stat(3,1) = rZlD;
fish_stat(3,2) = rmseZlD;
fish_stat(4,1) = rZlB;
fish_stat(4,2) = rmseZlB;
fish_stat(5,1) = rZ;
fish_stat(5,2) = rmseZ;
fish_stat(6,1) = rZl;
fish_stat(6,2) = rmseZl;
fish_stat(7,1) = rD;
fish_stat(7,2) = rmseD;
fish_stat(8,1) = rB;
fish_stat(8,2) = rmseB;

Fstat = array2table(fish_stat,'VariableNames',{'r','RMSE'},...
    'RowNames',{'ZD','ZB','ZlD','ZlB','Z','Zl','D','B'});
writetable(Fstat,[dpath 'LME_DvD_ZB_stats_' cfile '.csv'],'Delimiter',',','WriteRowNames',true)
save([dpath 'LME_DvD_ZB_stats_' cfile '.mat'],'fish_stat')

%% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

%% Figures
% Correlation
figure(1)
subplot(2,2,1)
for i=1:length(keep)
    lme=keep(i);
    plot(log10(RatPS(lme)),log10(RatZDet(lme)),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-0.25,-0.25,['r = ' num2str(rZD)])
xlabel('Fphotic:Fseabed')
ylabel('Zoop:Det')

subplot(2,2,2)
for i=1:length(keep)
    lme=keep(i);
    plot(log10(RatPS(lme)),log10(RatZB(lme)),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-0.25,-6.3,['r = ' num2str(rZB)])
xlabel('Fphotic:Fseabed')
ylabel('Zoop:Bent')

subplot(2,2,3)
for i=1:length(keep)
    lme=keep(i);
    plot(log10(RatPS(lme)),log10(RatZlDet(lme)),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-0.25,0.8,['r = ' num2str(rZlD)])
xlabel('Fphotic:Fseabed')
ylabel('ZoopLoss:Det')

subplot(2,2,4)
for i=1:length(keep)
    lme=keep(i);
    plot(log10(RatPS(lme)),log10(RatZlB(lme)),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-0.25,-5.2,['r = ' num2str(rZlB)])
xlabel('Fphotic:Fseabed')
ylabel('ZoopLoss:Bent')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_DvD_comp_ZBratios.png'])

%%
figure(2)
subplot(2,2,1)
for i=1:length(keep)
    lme=keep(i);
    plot(log10(Fphotic(lme)),log10(lme_z(lme)),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-0.45,3.5,['r = ' num2str(rZ)])
xlabel('Fphotic')
ylabel('Zoop')

subplot(2,2,3)
for i=1:length(keep)
    lme=keep(i);
    plot(log10(Fphotic(lme)),log10(lme_zl(lme)),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-0.45,4.75,['r = ' num2str(rZl)])
xlabel('Fphotic')
ylabel('ZoopLoss')

subplot(2,2,2)
for i=1:length(keep)
    lme=keep(i);
    plot(log10(Fseabed(lme)),log10(lme_det(lme)),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-1.75,4.75,['r = ' num2str(rD)])
xlabel('Fseabed')
ylabel('Det')

subplot(2,2,4)
for i=1:length(keep)
    lme=keep(i);
    plot(log10(Fseabed(lme)),log10(lme_mbio(lme,9)),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(-1.75,10.3,['r = ' num2str(rB)])
xlabel('Fseabed')
ylabel('Bent')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_DvD_comp_resources.png'])

%%
figure(3)
for i=1:length(keep)
    lme=keep(i);
    plot(log10(RatPS(lme)),log10(RatZlDet(lme)),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(log10(RatPS(lme)),log10(RatZlDet(lme)),num2str(lme),...
        'Color','k','HorizontalAlignment','center'); hold on;
end
text(-0.1,0.7,['r = ' num2str(rZlD)])
xlabel('Fphotic:Fseabed')
ylabel('ZoopLoss:Det')
axis([-0.2 1.8 -0.8 0.8])
stamp('')
print('-dpng',[ppath 'Clim_',harv,'_DvD_comp_ZlDet.png'])

figure(4)
for i=1:length(keep)
    lme=keep(i);
    plot(log10(RatPS(lme)),log10(RatZlB(lme)),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(log10(RatPS(lme)),log10(RatZlB(lme)),num2str(lme),...
        'Color','k','HorizontalAlignment','center'); hold on;
end
text(-0.1,-5.1,['r = ' num2str(rZlB)])
xlabel('Fphotic:Fseabed')
ylabel('ZoopLoss:Bent')
axis([-0.2 1.8 -6.5 -5])
stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_DvD_comp_ZlB.png'])


