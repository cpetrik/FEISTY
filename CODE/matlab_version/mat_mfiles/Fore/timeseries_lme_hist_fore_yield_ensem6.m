% Visualize time series output of FEISTY ensembles
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Saved as mat files
% Ensemble mid6, temp3
% Fishing yield in LMEs only
% total over all months in each year

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

%% Original parameters
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
fpath2 = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

%% Ensemble parameter sets w/baseline
efile = 'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050';

epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',efile,'/'];
epath2 = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/',efile,'/'];

load([epath 'LME_ts_yield_Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'Hlme_tsc_mf','Hlme_tsc_mp','Hlme_tsc_md','Hlme_tsc_lp','Hlme_tsc_ld');
load([epath 'LME_ts_yield_Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'Flme_tsc_mf','Flme_tsc_mp','Flme_tsc_md','Flme_tsc_lp','Flme_tsc_ld');

% Time series of catch (g per km2 per year)

%% ts
%In original saved file
y1 = 1861:2005;
y2 = 2006:2100;
y = [y1 y2];

HF = Hlme_tsc_mf;
HP = Hlme_tsc_mp + Hlme_tsc_lp;
HD = Hlme_tsc_md + Hlme_tsc_ld;
HLP = Hlme_tsc_lp;
HLD = Hlme_tsc_ld;

%% SOMETHING WRONG WITH MEAN OF #28 AROUND 1950
% FIX TS AT T=90 & 91 w/interp1
yi = y(89:92);
HF(28,89:92,61) = interp1(y([89,92]),HF(28,[89,92],61),yi);
HP(28,89:92,61) = interp1(y([89,92]),HP(28,[89,92],61),yi);
HD(28,89:92,61) = interp1(y([89,92]),HD(28,[89,92],61),yi);
HLP(28,89:92,61) = interp1(y([89,92]),HLP(28,[89,92],61),yi);
HLD(28,89:92,61) = interp1(y([89,92]),HLD(28,[89,92],61),yi);

%%
HA = HF + HP + HD;
HAA = HF + HLP + HLD;

FF = Flme_tsc_mf;
FP = Flme_tsc_mp + Flme_tsc_lp;
FD = Flme_tsc_md + Flme_tsc_ld;
FA = FF + FP + FD;
FLP = Flme_tsc_lp;
FLD = Flme_tsc_ld;
FAA = FF + FLP + FLD;

tF = [HF FF];
tP = [HP FP];
tD = [HD FD];
tA = [HA FA];
tLP = [HLP FLP];
tLD = [HLD FLD];
tAA = [HAA FAA];

%% difference from 1951
test=find(y>1950);
yid=test(1);

dtF = tF - tF(:,yid);
dtP = tP - tP(:,yid);
dtD = tD - tD(:,yid);
dtA = tA - tA(:,yid);
dtLP = tLP - tLP(:,yid);
dtLD = tLD - tLD(:,yid);
dtAA = tAA - tAA(:,yid);

% percent difference from 1951
pdF = (tF - tF(:,yid)) ./ tF(:,yid);
pdP = (tP - tP(:,yid)) ./ tP(:,yid);
pdD = (tD - tD(:,yid)) ./ tD(:,yid);
pdA = (tA - tA(:,yid)) ./ tA(:,yid);
pdLP = (tLP - tLP(:,yid)) ./ tLP(:,yid);
pdLD = (tLD - tLD(:,yid)) ./ tLD(:,yid);
pdAA = (tAA - tAA(:,yid)) ./ tAA(:,yid);

%% Cone of uncert raw diff
rmF = squeeze(mean(dtF));
rmP = squeeze(mean(dtP));
rmD = squeeze(mean(dtD));
rmA = squeeze(mean(dtA));
rmLP = squeeze(mean(dtLP));
rmLD = squeeze(mean(dtLD));
rmAA = squeeze(mean(dtAA));

rsF = squeeze(std(dtF));
rsP = squeeze(std(dtP));
rsD = squeeze(std(dtD));
rsA = squeeze(std(dtA));
rsLP = squeeze(std(dtLP));
rsLD = squeeze(std(dtLD));
rsAA = squeeze(std(dtAA));

%create continuous x value array for plotting
X=[y fliplr(y)];
%create y values for out and then back
%+/- 1 stdev
Ra=[rmA+rsA; flipud(rmA-rsA)];
Rf=[rmF+rsF; flipud(rmF-rsF)];
Rp=[rmP+rsP; flipud(rmP-rsP)];
Rd=[rmD+rsD; flipud(rmD-rsD)];
Raa=[rmAA+rsAA; flipud(rmAA-rsAA)];
Rlp=[rmLP+rsLP; flipud(rmLP-rsLP)];
Rld=[rmLD+rsLD; flipud(rmLD-rsLD)];

%% Cone pdiff
mpF = squeeze(mean(pdF));
mpP = squeeze(mean(pdP));
mpD = squeeze(mean(pdD));
mpA = squeeze(mean(pdA));
mpLP = squeeze(mean(pdLP));
mpLD = squeeze(mean(pdLD));
mpAA = squeeze(mean(pdAA));

spF = squeeze(std(pdF));
spP = squeeze(std(pdP));
spD = squeeze(std(pdD));
spA = squeeze(std(pdA));
spLP = squeeze(std(pdLP));
spLD = squeeze(std(pdLD));
spAA = squeeze(std(pdAA));

%create y values for out and then back
%+/- 1 stdev
Va=[mpA+spA; flipud(mpA-spA)];
Vf=[mpF+spF; flipud(mpF-spF)];
Vp=[mpP+spP; flipud(mpP-spP)];
Vd=[mpD+spD; flipud(mpD-spD)];
Vaa=[mpAA+spAA; flipud(mpAA-spAA)];
Vlp=[mpLP+spLP; flipud(mpLP-spLP)];
Vld=[mpLD+spLD; flipud(mpLD-spLD)];

%% Line color order
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

set(groot,'defaultAxesColorOrder',cm21);
% set(groot,'defaultFontName','TimesNewRoman');
% set(groot,'defaultFontSize',16);

nrows=11;
ncols=6;
pos = subfigrid(nrows,ncols,[0.05 0.025 0.025 0.05],[0.72 0.75]);

%scrsz = get(0,'ScreenSize');
% figure('units','pixels','position',...
%     [scrsz(3)/4,scrsz(4)/4,scrsz(3)/2,scrsz(4)/2])

%% all - baseline
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(y,(tA(end,:,i)*1e-6),'k','LineWidth',2); hold on;
        xlim([y(yid) y(end)])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Hist_Fore_',harv,'_yield_all_LMEs_subplot.png'])

%% types - baseline
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(y,(tF(end,:,i)*1e-6),'r','LineWidth',2); hold on;
        plot(y,(tP(end,:,i)*1e-6),'b','LineWidth',2); hold on;
        plot(y,(tD(end,:,i)*1e-6),'color',[0 0.6 0],'LineWidth',2); hold on;
        xlim([y(yid) y(end)])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end

print('-dpng',[ppath 'Hist_Fore_',harv,'_yield_types_LMEs_subplot.png'])

%% all baseline diff
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(y,(dtA(end,:,i)*1e-6),'k','LineWidth',2); hold on;
        xlim([y(yid) y(end)])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_yield_all_LMEs_subplot.png'])


%% types baseline diff
figure('Units','inches','Position',[1 3 6.5 8.5])
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(y,(dtF(end,:,i)*1e-6),'r','LineWidth',2); hold on;
        plot(y,(dtP(end,:,i)*1e-6),'b','LineWidth',2); hold on;
        plot(y,(dtD(end,:,i)*1e-6),'color',[0 0.6 0],'LineWidth',2); hold on;
        xlim([y(yid) y(end)])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_yield_types_LMEs_subplot.png'])

%% baseline adults diff
figure('Units','inches','Position',[1 3 6.5 8.5])
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(y,(dtAA(end,:,i)*1e-6),'k','LineWidth',2); hold on;
        xlim([y(yid) y(end)])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_yield_all_adults_LMEs_subplot.png'])

%% types - baseline diff
figure('Units','inches','Position',[1 3 6.5 8.5])
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(y,(dtF(end,:,i)*1e-6),'r','LineWidth',2); hold on;
        plot(y,(dtLP(end,:,i)*1e-6),'b','LineWidth',2); hold on;
        plot(y,(dtLD(end,:,i)*1e-6),'color',[0 0.6 0],'LineWidth',2); hold on;
        xlim([y(yid) y(end)])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_yield_types_adults_LMEs_subplot.png'])

%% adults w/uncert - diff
figure('Units','inches','Position',[1 3 6.5 8.5])
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        f=fill(X,(Raa(:,i))*1e-6,'k','FaceAlpha',0.25,'EdgeAlpha',0.25);  %plot filled area
        hold on
        plot(y,(rmAA(:,i))*1e-6,'k','LineWidth',2);
        xlim([y(yid) y(end)])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_yield_all_adults_LMEs_ensem_mid6_temp3_cone_1std_yr.png'])

%% adult types w/uncert - diff
figure('Units','inches','Position',[1 3 6.5 8.5])
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        fill(X,(Rf(:,i))*1e-6,'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
        fill(X,(Rlp(:,i))*1e-6,'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
        fill(X,(Rld(:,i))*1e-6,'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
        plot(y,(rmF(:,i))*1e-6,'r','LineWidth',2); hold on;
        plot(y,(rmLP(:,i))*1e-6,'b','LineWidth',2); hold on;
        plot(y,(rmLD(:,i))*1e-6,'color',[0 0.6 0],'LineWidth',2); hold on;
        xlim([y(yid) y(end)])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_yield_types_adults_LMEs_ensem_mid6_temp3_cone_1std_yr.png'])




