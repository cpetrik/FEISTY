% Time series of TEeff
% Need: AllL mdet  mmz_loss  mlz_loss mnpp

clear all
close all

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
gpath='/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

% Zoop and det and npp 
load([gpath 'cobalt_det_temp_zoop_npp_ts.mat']);

% FEISTY
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
BE = 0.075;
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'Means_fore_',harv,'_' cfile '.mat'],'ForeFish');
load([fpath 'Means_Historic_',harv,'_' cfile '.mat'],'HistFish');

cmYOR=cbrewer('seq','YlOrRd',50);
cmRP=cbrewer('seq','RdPu',50);
cmPR=cbrewer('seq','PuRd',50);


%% Zoop and det and npp 

%ESM2M in mmol N m-2 or mmol N m-2 d-1
% molN/m2 --> g/m2
% 106/16 mol C in 1 mol N
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W
mz_mean_fore = mz_mean_fore * (106.0/16.0) * 12.01 * 9.0;
lz_mean_fore = lz_mean_fore * (106.0/16.0) * 12.01 * 9.0;
% molN/m2/s --> g/m2/d
mzloss_mean_fore = mzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_mean_fore = lzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_mean_fore = det_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_mean_fore = npp_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;


mmz_mean = mz_mean_fore;
mlz_mean = lz_mean_fore;
mmz_loss = mzloss_mean_fore;
mlz_loss = lzloss_mean_fore;
mdet = det_mean_fore;
mnpp = npp_mean_fore;

% tmz_mean = mz_tot_hist * 1e-3 * 9.0;
% tlz_mean = lz_tot_hist * 1e-3 * 9.0;
% tmz_loss = mzl_tot_hist * 1e-3 * 9.0;
% tlz_loss = lzl_tot_hist * 1e-3 * 9.0;
% tdet = det_tot_hist * 1e-3 * 9.0;
% tnpp = npp_tot_hist * 1e-3 * 9.0;

%% plot info

geolon_t=double(geolon_t);
geolat_t=double(geolat_t);
[ni,nj]=size(geolon_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

Psf=NaN*ones(ni,nj);
Psp=NaN*ones(ni,nj);
Psd=NaN*ones(ni,nj);
Pmf=NaN*ones(ni,nj);
Pmp=NaN*ones(ni,nj);
Pmd=NaN*ones(ni,nj);
Plp=NaN*ones(ni,nj);
Pld=NaN*ones(ni,nj);
Pb=NaN*ones(ni,nj);

Psf(ID)=sf_prod50;
Psp(ID)=sp_prod50;
Psd(ID)=sd_prod50;
Pmf(ID)=mf_prod50;
Pmp(ID)=mp_prod50;
Pmd(ID)=md_prod50;
Plp(ID)=lp_prod50;
Pld(ID)=ld_prod50;
Pb(ID)=b_mean50;

Psf(Psf(:)<0) = 0;
Psp(Psp(:)<0) = 0;
Psd(Psd(:)<0) = 0;
Pmf(Pmf(:)<0) = 0;
Pmp(Pmp(:)<0) = 0;
Pmd(Pmd(:)<0) = 0;
Plp(Plp(:)<0) = 0;
Pld(Pld(:)<0) = 0;
Pb(Pb(:)<0) = 0;

All = Psp+Psf+Psd+Pmp+Pmf+Pmd+Plp+Pld;
AllF = Psf+Pmf;
AllP = Psp+Pmp+Plp;
AllD = Psd+Pmd+Pld;
AllS = Psp+Psf+Psd;
AllM = Pmp+Pmf+Pmd;
AllL = Plp+Pld;

%%
figure
subplot(2,2,1)
hist(log10(mnpp(:)))
title('NPP')
subplot(2,2,2)
hist(log10(mdet(:)))
title('Det')
subplot(2,2,3)
hist(log10(mmz_loss(:) + mlz_loss(:)))
title('Zoop')
subplot(2,2,4)
hist(log10(AllL(:)))
title('Large')

%% Effective TEs
% With BE*det instead of Bent
TEeffM = AllM./(BE*mdet + mmz_loss + mlz_loss); 
%TEeff_L = production_L/NPP
TEeff_L = AllL./mnpp;
TEeff_L(TEeff_L==-Inf) = NaN;
TEeff_L(TEeff_L==Inf) = NaN;
TEeff_L(TEeff_L<0) = NaN;
%TEeff_LTL = (production_benthic_invert+mesozoo_prod_to_fish)/NPP
TEeff_LTLd = (BE*mdet + mmz_loss + mlz_loss)./mnpp;
TEeff_LTLd(TEeff_LTLd==-Inf) = NaN;
TEeff_LTLd(TEeff_LTLd==Inf) = NaN;
TEeff_LTLd(TEeff_LTLd<0) = NaN;
%TEeff_HTL = production_L/(production_benthic_invert+mesozoo_prod_to_fish)
TEeff_HTLd = AllL./(BE*mdet + mmz_loss + mlz_loss); 
TEeff_HTLd(TEeff_HTLd<0) = NaN;

TELTLd1 = real(TEeff_LTLd.^(1/1.25));
TELTLd2 = real(TEeff_LTLd.^(1/1.5));

TEM = real(TEeffM.^(1/2));          %should this be 1/1?
TEL = real(TEeff_L.^(1/4));         %should this be 1/3?
TEHTLd = real(TEeff_HTLd.^(1/3));   %should this be 1/2?

q(:,1) = [0.01 0.05 0.25 0.5 0.75 0.95 0.99]';
q(:,2) = quantile((TEeff_LTLd(:)),[0.01 0.05 0.25 0.5 0.75 0.95 0.99])';
q(:,3) = quantile((TEeff_HTLd(:)),[0.01 0.05 0.25 0.5 0.75 0.95 0.99]);
q(:,4) = quantile((TEeff_L(:)),[0.01 0.05 0.25 0.5 0.75 0.95 0.99]);
q(:,5) = quantile((TEHTLd(:)),[0.01 0.05 0.25 0.5 0.75 0.95 0.99]);
q(:,6) = quantile((TEL(:)),[0.01 0.05 0.25 0.5 0.75 0.95 0.99]);


Q = array2table(q,'VariableNames',{'Quantile','TEeff_LTLd','TEeff_HTLd',...
    'TEeff_L','TEHTLd','TEL'});

%% save
writetable(Q,[fpath 'TEeff_quant_Forecast_All_fish03_' cfile '.csv'],'Delimiter',',');

save([fpath 'TEeffDet_Forecast_All_fish03_' cfile '.mat'],'TEeffM',...
    'Pmf','Pmp','Pmd','Plp','Pld','Pb','mmz_loss','mlz_loss','mnpp',...
    'TEeff_L','TEeff_LTLd','TEeff_HTLd');
