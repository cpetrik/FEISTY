% Compare the two grids to see if they are different

clear all
close all

fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/advect_tests/';

load('/Volumes/GFDL/GCM_DATA/CORE-forced/feb152013_run25_ocean.198801-200712_uh200_vh200.mat',...
    'uh200','vh200','u200','v200');

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'AREA_OCN','dat','dxtn','dyte','ht','geolon_t','geolat_t');


area = AREA_OCN;
area = (area)*510072000*1e6;
area = max(area,1);

% grid size
ni = 360;
nj = 200;
isd = 1;
jsd = 1;
ied = ni;
jed = nj;

% depth of the surface layer, 200m or less
eps = 1;
dep = min(ht,200);
dep = max(dep,eps);
% dep = ones(ni,nj);

% land mask
mask = zeros(ni,nj);
aa = find(ht > 0);
mask(aa) = 1;

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Data_hindcast_grid_cp2D.mat')

dxtn2 = GRD.dxtn;
dyte2 = GRD.dyte;
area2 = GRD.dat;
mask2 = GRD.mask;
dep2 = GRD.ht;

area2 = area2.*mask2;
area2 = max(area2,1);

dep2 = min(dep2,200);
dep2 = max(dep2,1);

dx = (dxtn~=dxtn2);
dy = (dyte~=dyte2);
ar = (area~=area2);
ma = (mask~=mask2);
dp = (dep~=dep2);

sum(dx(:))
sum(dy(:))
sum(ar(:))
sum(ma(:))
sum(dp(:))






