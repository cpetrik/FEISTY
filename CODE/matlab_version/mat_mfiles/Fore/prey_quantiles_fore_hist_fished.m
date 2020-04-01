% Calc prey quantiles (for scope for growth calc)
% ESM2M Hindcast of 1951-2000 and Forecast of 2051-2100

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% NPP and zoop
load([bpath 'cobalt_det_temp_zoop_npp_means.mat']);

% molN/m2/s --> g/m2/d
mz_hist = mz_mean_hist * (106.0/16.0) * 12.01 * 9.0;
lz_hist = lz_mean_hist * (106.0/16.0) * 12.01 * 9.0;

mz_fore = mz_mean_fore * (106.0/16.0) * 12.01 * 9.0;
lz_fore = lz_mean_fore * (106.0/16.0) * 12.01 * 9.0;

z_hist = mz_hist + lz_hist;
z_fore = mz_fore + lz_fore;


%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
%fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

% Hindcast
load([fpath 'Means_Historic_' harv '_' cfile '.mat'],...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50','b_mean50');

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
Hsf(grid(:,1))=sf_mean50;
Hsp(grid(:,1))=sp_mean50;
Hsd(grid(:,1))=sd_mean50;
Hmf(grid(:,1))=mf_mean50;
Hmp(grid(:,1))=mp_mean50;
Hmd(grid(:,1))=md_mean50;
Hlp(grid(:,1))=lp_mean50;
Hld(grid(:,1))=ld_mean50;
Hb(grid(:,1)) =b_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50 b_mean50

% Forecast
load([fpath 'Means_fore_',harv,'_' cfile '.mat'],...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50','b_mean50');

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
Csf(ID)=sf_mean50;
Csp(ID)=sp_mean50;
Csd(ID)=sd_mean50;
Cmf(ID)=mf_mean50;
Cmp(ID)=mp_mean50;
Cmd(ID)=md_mean50;
Clp(ID)=lp_mean50;
Cld(ID)=ld_mean50;
Cb(ID) =b_mean50;

clear sf_mean50 sp_mean50 sd_mean50 mf_mean50 mp_mean50 md_mean50 lp_mean50 ld_mean50 b_mean50

cF = Csf+Cmf;
cP = Csp+Cmp+Clp;
cD = Csd+Cmd+Cld;
cS = Csp+Csf+Csd;
cM = Cmp+Cmf+Cmd;
cL = Clp+Cld;

hF = Hsf+Hmf;
hP = Hsp+Hmp+Hlp;
hD = Hsd+Hmd+Hld;
hS = Hsp+Hsf+Hsd;
hM = Hmp+Hmf+Hmd;
hL = Hlp+Hld;

cAll = cF+cP+cD;
hAll = hF+hP+hD;

%% plot info
plotminlat=60; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%% Prey
% Sm = 0.25;  %Feeding 2 sizes down
% J = 1.0;    %Juvenile feeding reduction
% D = 0.75;   %Demersal feeding in pelagic reduction
% A = 0.5;    %Adult predation reduction %*****
% 
% MF_phi_MZ = Sm;
% MF_phi_LZ = 1.0;
% MF_phi_S  = 1.0;
% 
% MP_phi_MZ = Sm*J;
% MP_phi_LZ = J;
% MP_phi_S  = J;
% 
% MD_phi_BE = 1.0;
% 
% LP_phi_MF = 1.0*A;
% LP_phi_MP = 1.0;
% LP_phi_MD = 0.0;
% 
% LD_phi_MF = D*A;
% LD_phi_MP = D;
% LD_phi_MD = 1.0;
% LD_phi_BE = 1.0;

Mpel_hist = hS + lz_hist + 0.25*mz_hist;
Mpel_fore = cS + lz_fore + 0.25*mz_fore;

Mdem_hist = Hb;
Mdem_fore = Cb;

Lpel_hist = hM;
Lpel_fore = cM;

Ldem_hist = Hb + 0.75*hM;
Ldem_fore = Cb + 0.75*cM;

M_all = [Mpel_hist Mpel_fore Mdem_hist Mdem_fore];
L_all = [Lpel_hist Lpel_fore Ldem_hist Ldem_fore];
all = [M_all L_all];

%% quantiles
q(:,1) = [0.05 0.25 0.5 0.75 0.95]';
q(:,2) = quantile((Mpel_hist(:)),[0.05 0.25 0.5 0.75 0.95]);
q(:,3) = quantile((Mpel_fore(:)),[0.05 0.25 0.5 0.75 0.95]);
q(:,4) = quantile((Mdem_hist(:)),[0.05 0.25 0.5 0.75 0.95]);
q(:,5) = quantile((Mdem_fore(:)),[0.05 0.25 0.5 0.75 0.95]);
q(:,6) = quantile((Lpel_hist(:)),[0.05 0.25 0.5 0.75 0.95]);
q(:,7) = quantile((Lpel_fore(:)),[0.05 0.25 0.5 0.75 0.95]);
q(:,8) = quantile((Ldem_hist(:)),[0.05 0.25 0.5 0.75 0.95]);
q(:,9) = quantile((Ldem_fore(:)),[0.05 0.25 0.5 0.75 0.95]);
q(:,10) = quantile((M_all(:)),[0.05 0.25 0.5 0.75 0.95]);
q(:,11) = quantile((L_all(:)),[0.05 0.25 0.5 0.75 0.95]);
q(:,12) = quantile(all(:),[0.05 0.25 0.5 0.75 0.95]);

Q = array2table(q,'VariableNames',{'Quantile','Mpel_Hist','Mpel_Fore',...
    'Mdem_Hist','Mdem_Fore','Lpel_Hist','Lpel_Fore','Ldem_Hist','Ldem_Fore',...
    'M_all','L_all','All'});

writetable(Q,[fpath 'Prey_quant_Fore_Hist_All_fish03_' cfile '.csv'],'Delimiter',',');

save([fpath 'Prey_quant_Fore_Hist_All_fish03_' cfile '.mat'],'q','Q');



