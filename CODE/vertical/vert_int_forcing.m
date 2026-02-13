% Vertically integrated forcing for orig FEISTY

clear
close all

%%
vpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
%vpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%1-D
load([vpath 'grid_OM4_05_COBALTv3.mat'],'wet','geolon','z_l');

%ocean grid cells
WID = find(wet(:)==1);
NWID = length(WID);
NZID = length(z_l);
[ni,nj] = size(geolon);

%%
[TEMP_z,TEMP_btm,MZ,LZ,MZloss,LZloss,det_btm] = ncread_global_feisty_forcing_first_year(vpath,ni,nj,NZID);

% THICKNESS
load([vpath 'ocean_cobalt_feisty_forcing_z.199001-199412.thkcello.mat'],'thkcello')
thkcello = thkcello(:,:,:,1:12);

%%
Time=15:30:365;
Tdays=1:365;

%%

%Interpolate monthly forcing to daily
ESM = daily_interp_int_monthly_means(NWID,Time,Tdays,...
    TEMP_z,TEMP_btm,det_btm,MZ,LZ,MZloss,LZloss,thkcello,WID,ni,nj);

%% Save
save([vpath 'ocean_cobalt_feisty_forcing_2dint.1990.mat'],'ESM','-v7.3')