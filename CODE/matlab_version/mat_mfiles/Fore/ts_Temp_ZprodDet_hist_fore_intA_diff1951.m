% Time series of P vs D and Z vs Det
% Area-integrated totals
% Z prod instead off loss

clear 
close all

%%
cm={[0.5 0 1],...   %purple
    [1 0 1],...     %m
    [1 0 0],...     %r
    [0.5 0 0],...   %maroon
    [0.75 0.75 0.75],... %lt grey
    [0.5 0.5 0.5],...    %med grey
    [49/255 79/255 79/255],... %dk grey
    [0 0 0],...      %black
    [1 1 0],...      %yellow
    [127/255 255/255 0],... %lime green
    [0 0.5 0],...    %dk green
    [0/255 206/255 209/255],... %turq
    [0 0.5 0.75],...   %med blue
    [188/255 143/255 143/255],... %rosy brown
    [255/255 192/255 203/255],... %pink
    [255/255 160/255 122/255]}; %peach

%% ts
y1 = 1860+(1/12):5:2005;
y2 = 2005+(1/12):5:2100;
y = [y1 y2];

%% Zoop:Det
% Visualize difference between ESM2M Hindcast (1951-2000) and Forecast (2051-2100)

cpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/cobalt_data/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ffold=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/'];
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

%% Zoop & Det
load([bpath 'cobalt_det_temp_zoop_npp_means.mat'],'mzprod_1yr_hist',...
    'mzprod_1yr_fore','lzprod_1yr_hist','lzprod_1yr_fore','det_1yr_hist',...
    'det_1yr_fore','ptemp_1yr_hist','ptemp_1yr_fore',...
    'ptemp_5yr_hist','ptemp_5yr_fore');
load([ffold 'ESM2M_Hist_Fore/ts_Hist_Fore_Zp_D_ZpDet_intA.mat']);

% molN/m2/s --> g/m2/d
mzprod_hist = mzprod_1yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_hist = lzprod_1yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_hist = det_1yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

mzprod_fore = mzprod_1yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_fore = lzprod_1yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_fore = det_1yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

zprod_hist = mzprod_hist + lzprod_hist;
zprod_fore = mzprod_fore + lzprod_fore;


%%
y1 = 1860+(1/12):1:2005;
y2 = 2005+(1/12):1:2100;
yr = [y1 y2];
id = find(yr>1950);
yid = id(1);

% times series of 1-yr (not area averaged?)
ptemp1 = [ptemp_1yr_hist ptemp_1yr_fore] - 273;
ptemp5 = [ptemp_5yr_hist ptemp_5yr_fore] - 273;

det = [det_hist det_fore];
zprod = [zprod_hist zprod_fore];

mtp1 = mean(ptemp1);
mtp5 = mean(ptemp5);
mdet = mean(det);
mzp = mean(zprod);

dT1 = mtp1 - mtp1(yid);

pdiffDet = (mdet-mdet(yid)) ./ mdet(yid);
pdiffZ = (mzp-mzp(yid)) ./ mzp(yid);

%time series of 5-yr 1yrs as difference from 1951
dtD = tD - tD(19);
dtZ = tZ - tZ(19);
dT5 = mtp5 - mtp5(19);

pdtD = (tD - tD(19)) / tD(19);
pdtZ = (tZ - tZ(19)) / tZ(19);


%% ts 
figure(1)
subplot(2,2,1)
line(y(19:end),dT5(19:end),'color','k','Linewidth',2); hold on;
ylabel('Change in temperature (^oC)');

subplot(2,2,3)
line(y(19:end),100*pdtD(19:end),'color',[0 0 0.45],'Linewidth',2); hold on;
line(y(19:end),100*pdtZ(19:end),'color',[0.35 0 0.6],'Linewidth',2); hold on;
ylabel('Percent change in production');
xlabel('Year')
print('-dpng',[pp 'Hist_Fore_Zp_Det_temp_diff1951.png'])

%% Zprod, Det not area weighted
figure (2)
subplot(2,2,1)
line(yr(yid:end),dT1(yid:end),'color','k','Linewidth',2); hold on;

subplot(2,2,3)
line(yr(yid:end),pdiffDet(yid:end),'color',cm{1},'Linewidth',2); hold on;
line(yr(yid:end),pdiffZ(yid:end),'color',cm{12},'Linewidth',2); hold on;
xlabel('Year');
%print('-dpng',[pp 'Hist_Fore_',harv,'_PDfrac_ZpDet_temp_diff1951.png'])

