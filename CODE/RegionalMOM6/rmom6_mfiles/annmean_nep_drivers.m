% View time series and maps
% Reanalysis-forced runs 1993-2023, v20250715
% MOM6-NEP10

clear
close all

cpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-NEP10/';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/',...
    'FEISTY/CODE/Figs/Matlab_New_sizes/',...
    'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/',...
    'NEP10/'];

%% Units
%poc flux: mol N m-2 s-1
%zoo: mol N m-2
%tp: degC
%tb: degC

load([cpath 'temp_100.nep.full.hcast.monthly.raw.r20250912.199301-202506.mat'],'temp_100');
load([cpath 'btm_temp.nep.full.hcast.monthly.raw.r20250912.199301-202506.mat'],'btm_temp');
load([cpath 'fntot_btm.nep.full.hcast.monthly.raw.r20250912.199301-202506.mat'],'fntot_btm','fntot_btm_long_name','fntot_btm_units')
load([cpath 'nmdz_nlgz_100.nep.full.hcast.monthly.raw.r20250912.199301-202506.mat']);
load([cpath 'jhploss_mdz_lgz_100.nep.full.hcast.monthly.raw.r20250912.199301-202506.mat']);

%% nans
temp_100(temp_100 > 1.0e19) = nan;
btm_temp(btm_temp > 1.0e19) = nan;
fntot_btm(fntot_btm > 1.0e19) = nan;
nmdz_100(nmdz_100 > 1.0e19) = nan;
nlgz_100(nlgz_100 > 1.0e19) = nan;
jhploss_nmdz_100(jhploss_nmdz_100 > 1.0e19) = nan;
jhploss_nlgz_100(jhploss_nlgz_100 > 1.0e19) = nan;

%% units
nmdz_100 = nmdz_100 * (106/16) * 12.01  * 9.0;
nlgz_100 = nlgz_100 * (106/16) * 12.01  * 9.0;
jhploss_nmdz_100 = jhploss_nmdz_100 * (106/16) * 12.01 * 9.0 * 60 * 60 * 24;
jhploss_nlgz_100 = jhploss_nlgz_100 * (106/16) * 12.01 * 9.0 * 60 * 60 * 24;
fntot_btm = fntot_btm * (106/16) * 12.01 * 9.0 * 60 * 60 * 24;

%%
load([cpath 'nep_raw_ocean_static_gridspec.mat'],'geolon','geolat');
load([cpath 'Data_grid_mom6_nep10.mat'],'GRD');

[ni,nj]=size(geolon);
geolon = double(geolon);
geolat = double(geolat);

%NE Pac
plotminlat=10; 
plotmaxlat=85;
plotminlon=156;
plotmaxlon=-104;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;   %decent looking coastlines

%% colors
cm10=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black
    

set(groot,'defaultAxesColorOrder',cm10);

%% time
%units = 'days since 1993-01-01 00:00:00'
mos = length(time);
nyrs = floor(mos/12);

%yrs = floor(yr(1)):floor(yr(end));
yrs = floor(yr(1)):(floor(yr(1))+nyrs-1);

%% annual means

st=1:12:mos;
en=12:12:mos;

tp_amean = nan*ones(ni,nj,nyrs);
tb_amean = tp_amean;
mz_amean = tp_amean;
lz_amean = tp_amean;
mzhp_amean = tp_amean;
lzhp_amean = tp_amean;
det_amean = tp_amean;

for n=1:nyrs
    % mean 
    tp_amean(:,:,n) = squeeze(mean(temp_100(:,:,st(n):en(n)),3,'omitnan'));
    tb_amean(:,:,n) = squeeze(mean(btm_temp(:,:,st(n):en(n)),3,'omitnan'));
    mz_amean(:,:,n) = squeeze(mean(nmdz_100(:,:,st(n):en(n)),3,'omitnan'));
    lz_amean(:,:,n) = squeeze(mean(nlgz_100(:,:,st(n):en(n)),3,'omitnan'));
    mzhp_amean(:,:,n) = squeeze(mean(jhploss_nmdz_100(:,:,st(n):en(n)),3,'omitnan'));
    lzhp_amean(:,:,n) = squeeze(mean(jhploss_nlgz_100(:,:,st(n):en(n)),3,'omitnan'));
    det_amean(:,:,n) = squeeze(mean(fntot_btm(:,:,st(n):en(n)),3,'omitnan'));

end

%% seasonal cycle
tp_smean = nan*ones(ni,nj,12);
tb_smean = tp_smean;
mz_smean = tp_smean;
lz_smean = tp_smean;
mzhp_smean = tp_smean;
lzhp_smean = tp_smean;
det_smean = tp_smean;

for n=1:12
    s = n:12:mos;
    % mean 
    tp_smean(:,:,n)=mean(temp_100(:,:,s),3,'omitnan');
    tb_smean(:,:,n)=mean(btm_temp(:,:,st(n):en(n)),3,'omitnan');
    mz_smean(:,:,n)=mean(nmdz_100(:,:,st(n):en(n)),3,'omitnan');
    lz_smean(:,:,n)=mean(nlgz_100(:,:,st(n):en(n)),3,'omitnan');
    mzhp_smean(:,:,n)=mean(jhploss_nmdz_100(:,:,st(n):en(n)),3,'omitnan');
    lzhp_smean(:,:,n)=mean(jhploss_nlgz_100(:,:,st(n):en(n)),3,'omitnan');
    det_smean(:,:,n)=mean(fntot_btm(:,:,st(n):en(n)),3,'omitnan');

end


%% SAVE
save([cpath 'nep.full.hcast.monthly.raw.r20250912.199301-202506_annual_means.mat'],'yrs',...
    'tp_amean','tb_amean','mz_amean','lz_amean','mzhp_amean','lzhp_amean','det_amean',...
    'tp_smean','tb_smean','mz_smean','lz_smean','mzhp_smean','lzhp_smean','det_smean')

