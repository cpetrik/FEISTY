% Bar plot of trophic amplification

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t','AREA_OCN');
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);
load([cpath 'lme_mask_esm2m.mat']);

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

%% NPP and zoop
load([bpath 'cobalt_det_temp_zoop_npp_means.mat']);

% molN/m2/s --> g/m2/d
mzloss_hist = mzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_hist = lzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_hist = npp_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_hist = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
ptemp_hist = ptemp_mean_hist - 273;

mzloss_fore = mzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_fore = lzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_fore = npp_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_fore = det_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
ptemp_fore = ptemp_mean_fore - 273;

zloss_hist = mzloss_hist + lzloss_hist;
zloss_fore = mzloss_fore + lzloss_fore;

% g/m2/d --> total g; Mult flux by 365 - YES
Azloss_hist = 365 * zloss_hist .* AREA_OCN;
Azloss_fore = 365 * zloss_fore .* AREA_OCN;
Adet_hist = 365 * det_hist .* AREA_OCN;
Adet_fore = 365 * det_fore .* AREA_OCN;
Anpp_hist = 365 * npp_hist .* AREA_OCN; 
Anpp_fore = 365 * npp_fore .* AREA_OCN; 
%nansum(Anpp_mean(:))/9 ~ 50 PgC (1PgC = 1e15gC)
%nansum(Anpp_hist(:))/9 = 5.4743e+16

%grams to tons
%1e-6*5e16 = 5.0000e+10
%1.10231e-6*5e16 = 5.5116e+10

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

% Hindcast
load([fpath 'Means_Historic_' harv '_prod_' cfile '.mat'],...
    'sf_prod50','sp_prod50','sd_prod50',...
    'mf_prod50','mp_prod50','md_prod50',...
    'lp_prod50','ld_prod50','b_mean50');

[hi,hj]=size(geolon_t);
Hsf=NaN*ones(hi,hj);
Hsp=NaN*ones(hi,hj);
Hsd=NaN*ones(hi,hj);
Hmf=NaN*ones(hi,hj);
Hmp=NaN*ones(hi,hj);
Hmd=NaN*ones(hi,hj);
Hlp=NaN*ones(hi,hj);
Hld=NaN*ones(hi,hj);
Hb =NaN*ones(hi,hj);

Hsf(grid(:,1))=sf_prod50;
Hsp(grid(:,1))=sp_prod50;
Hsd(grid(:,1))=sd_prod50;
Hmf(grid(:,1))=mf_prod50;
Hmp(grid(:,1))=mp_prod50;
Hmd(grid(:,1))=md_prod50;
Hlp(grid(:,1))=lp_prod50;
Hld(grid(:,1))=ld_prod50;
Hb(grid(:,1)) =b_mean50;

Hsf(Hsf(:)<0)=0;
Hsp(Hsp(:)<0)=0;
Hsd(Hsd(:)<0)=0;
Hmf(Hmf(:)<0)=0;
Hmp(Hmp(:)<0)=0;
Hmd(Hmd(:)<0)=0;
Hlp(Hlp(:)<0)=0;
Hld(Hld(:)<0)=0;

Hsf_area = Hsf .* AREA_OCN;
Hsp_area = Hsp .* AREA_OCN;
Hsd_area = Hsd .* AREA_OCN;
Hmf_area = Hmf .* AREA_OCN;
Hmp_area = Hmp .* AREA_OCN;
Hmd_area = Hmd .* AREA_OCN;
Hlp_area = Hlp .* AREA_OCN;
Hld_area = Hld .* AREA_OCN;
Hb_area = Hb .* AREA_OCN;

clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50

%% Forecast
load([fpath 'Means_fore_',harv,'_' cfile '.mat'],...
    'sf_prod50','sp_prod50','sd_prod50',...
    'mf_prod50','mp_prod50','md_prod50',...
    'lp_prod50','ld_prod50','b_mean50');

[ni,nj]=size(geolon_t);
Csf=NaN*ones(ni,nj);
Csp=NaN*ones(ni,nj);
Csd=NaN*ones(ni,nj);
Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);
Cb =NaN*ones(ni,nj);

Csf(ID)=sf_prod50;
Csp(ID)=sp_prod50;
Csd(ID)=sd_prod50;
Cmf(ID)=mf_prod50;
Cmp(ID)=mp_prod50;
Cmd(ID)=md_prod50;
Clp(ID)=lp_prod50;
Cld(ID)=ld_prod50;
Cb(ID) =b_mean50;

Csf(Csf(:)<0)=0;
Csp(Csp(:)<0)=0;
Csd(Csd(:)<0)=0;
Cmf(Cmf(:)<0)=0;
Cmp(Cmp(:)<0)=0;
Cmd(Cmd(:)<0)=0;
Clp(Clp(:)<0)=0;
Cld(Cld(:)<0)=0;

Csf_area = Csf .* AREA_OCN;
Csp_area = Csp .* AREA_OCN;
Csd_area = Csd .* AREA_OCN;
Cmf_area = Cmf .* AREA_OCN;
Cmp_area = Cmp .* AREA_OCN;
Cmd_area = Cmd .* AREA_OCN;
Clp_area = Clp .* AREA_OCN;
Cld_area = Cld .* AREA_OCN;
Cb_area = Cb .* AREA_OCN;

clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50

%% Groups
cF = Csf_area + Cmf_area;
cP = Csp_area + Cmp_area + Clp_area;
cD = Csd_area + Cmd_area + Cld_area;
cS = Csp_area + Csf_area + Csd_area;
cM = Cmp_area + Cmf_area + Cmd_area;
cL = Clp_area + Cld_area;
cAll = cF+cP+cD;
cPel = cF+cP;

hF = Hsf_area + Hmf_area;
hP = Hsp_area + Hmp_area + Hlp_area ;
hD = Hsd_area + Hmd_area + Hld_area ;
hS = Hsp_area + Hsf_area + Hsd_area ;
hM = Hmp_area + Hmf_area + Hmd_area ;
hL = Hlp_area + Hld_area;
hAll = hF+hP+hD;
hPel = hF+hP;

%% Percent differences

pbar(1) = (nansum(Anpp_fore(:))-nansum(Anpp_hist(:))) ./ nansum(Anpp_hist(:));
pbar(2) = -0.079; %(nansum(Azloss_fore(:))-nansum(Azloss_hist(:))) ./ nansum(Azloss_hist(:));
pbar(3) = (nansum(Adet_fore(:))-nansum(Adet_hist(:))) ./ nansum(Adet_hist(:));
pbar(4) = (nansum(cM(:))-nansum(hM(:))) ./ nansum(hM(:));
pbar(5) = (nansum(cL(:))-nansum(hL(:))) ./ nansum(hL(:));
pbar(6) = (nansum(cF(:))-nansum(hF(:))) ./ nansum(hF(:));
pbar(7) = (nansum(cP(:))-nansum(hP(:))) ./ nansum(hP(:));
pbar(8) = (nansum(cD(:))-nansum(hD(:))) ./ nansum(hD(:));
pbar(9) = (nansum(cPel(:))-nansum(hPel(:))) ./ nansum(hPel(:));
pbar(10) = (nansum(cAll(:))-nansum(hAll(:))) ./ nansum(hAll(:));
pnames = {'NPP','Z','Det','M','L','F','P','D','Pel','All'};

delB = (nansum(Cb_area(:))-nansum(Hb_area(:))) ./ nansum(Hb_area(:)); 
%change in B is -10.16%, but change in Det is -12.3%

%% bar graphs
figure(1)
bar(100*pbar)
set(gca,'XTickLabel',pnames)

%% size
figure(2)
bar(pbar([1,2,4,5])*100)
ylim([-25 0])
colormap('gray')
set(gca,'XTickLabel',{'NPP','Mesozoo','M Fishes','L Fishes'})
ylabel('Percent change')
title('Global change in productivity')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_size.png'])

figure(3)
bar(pbar([1,2,4,5])*100)
ylim([-25 0])
colormap('gray')
set(gca,'XTickLabel',[])
% ylabel('Percent change')
% title('Global change in productivity')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_size_nolabs.png'])

%% type
figure(4)
bar(pbar([6:8])*100)
ylim([-25 0])
colormap('gray')
set(gca,'XTickLabel',{'F','P','D'})
ylabel('Percent change')
title('Global change in productivity')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_type.png'])



