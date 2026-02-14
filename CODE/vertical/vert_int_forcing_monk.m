% Vertically integrated forcing for orig FEISTY

clear
close all

%%
vpath = '/scratch/cpetrik/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
%vpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%1-D
load([vpath 'grid_OM4_05_COBALTv3.mat'],'wet','geolon','z_l');

%ocean grid cells
WID = find(wet(:)==1);
NWID = length(WID);
NZID = length(z_l);
[ni,nj] = size(geolon);

%% years
ystart = 1990:5:2015;
yend = 1994:5:2019;
nyears = length(ystart);

Time=15:30:365;
Tdays=1:365;

%% loop over years
for y = 1:nyears

    %%
    [TEMP_z,TEMP_btm,MZ,LZ,MZloss,LZloss,det_btm] = ncread_global_feisty_forcing_first_year(vpath,ni,nj,NZID);

    % THICKNESS
    load([vpath 'ocean_cobalt_feisty_forcing_z.',num2str(ystart(y)),'01-',num2str(yend(y)),'12.thkcello.mat'],'thkcello')
    thkcello = thkcello(:,:,:,1:12);

    %% Interpolate monthly forcing to daily
    ESM = daily_interp_int_monthly_means(NWID,Time,Tdays,...
        TEMP_z,TEMP_btm,det_btm,MZ,LZ,MZloss,LZloss,thkcello,WID,ni,nj);

    %% Save
    save([vpath 'ocean_cobalt_feisty_forcing_2dint_',num2str(ystart(y)),'.mat'],'ESM',"-v7.3")

end