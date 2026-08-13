% Make grid data for NEP10
% TP and vel have fewer ocean cells than BGC
% Comp to wet?

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-NEP10/';

%%
load([fpath 'nep_raw_ocean_static_gridspec.mat'])

%% index of water cells
WID = find((wet(:)==1));
WIDc = find((wet_c(:)==1));
WIDu = find((wet_u(:)==1));
WIDv = find((wet_v(:)==1));

WIDtu = intersect(WID,WIDu);
WIDtuv = intersect(WIDv,WIDtu);

% WID         205698x1
% WIDc        204122x1
% WIDtu       152495x1
% WIDtuv      151971x1
% WIDu        204994x1
% WIDv        204892x1

%%
NID = length(WID);

GRD.ID = WID;
GRD.NID = NID;
GRD.Lat = double(geolat(WID));
GRD.Lon = double(geolon(WID));
GRD.depth = double(deptho(WID));
GRD.area = double(areacello(WID));

save([fpath 'Data_grid_mom6_nep10.mat'], 'GRD');

%%
[ni,nj] = size(geolon);
mask = zeros(ni,nj);
mask(WID) = ones;

GRD2.area = double(areacello);
GRD2.dx = double(dxt);
GRD2.dy = double(dyt);
GRD2.lat = double(geolat);
GRD2.lon = double(geolon);
GRD2.mask = double(mask);

save([fpath 'Data_grid2D_mom6_nep10.mat'], 'GRD2');


