% Save PRODUCTION Hist and Fore together

clear all
close all

% Hindcast grid
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
lpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

harv = 'All_fish03';

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

clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50

% Forecast
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

clear sf_prod50 sp_prod50 sd_prod50 mf_prod50 mp_prod50 md_prod50 lp_prod50 ld_prod50 b_mean50

cFprod = Csf+Cmf;
cPprod = Csp+Cmp+Clp;
cDprod = Csd+Cmd+Cld;
cSprod = Csp+Csf+Csd;
cMprod = Cmp+Cmf+Cmd;
cLprod = Clp+Cld;

hFprod = Hsf+Hmf;
hPprod = Hsp+Hmp+Hlp;
hDprod = Hsd+Hmd+Hld;
hSprod = Hsp+Hsf+Hsd;
hMprod = Hmp+Hmf+Hmd;
hLprod = Hlp+Hld;

hB = Hb;
cB = Cb;

save([fpath 'Means_hist_fore_',harv,'_cobalt_' cfile '.mat'],'cFprod',...
    'cPprod','cDprod','cSprod','cMprod','cLprod','cB','hFprod',...
    'hPprod','hDprod','hSprod','hMprod','hLprod','hB','-append');
save([lpath 'Means_hist_fore_',harv,'_cobalt_' cfile '.mat'],'cFprod',...
    'cPprod','cDprod','cSprod','cMprod','cLprod','cB','hFprod',...
    'hPprod','hDprod','hSprod','hMprod','hLprod','hB','-append');

