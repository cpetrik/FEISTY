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
Time=15:30:365;
Tdays=1:365;

%% year indices for 5 yr, monthly means
st = 1:12:60;

yst = 1990:5:2019;
yen = 1994:5:2019;

for Y = 1:length(yst)
    %% read in one year
    yrs = [num2str(yst(Y)),'01-',num2str(yen(Y)),'12'];
    [TEMP_z,TEMP_btm,MZ,LZ,MZloss,LZloss,det_btm] = ncread_global_feisty_forcing_1year(vpath,ni,nj,NZID,yrs,styr);

    % THICKNESS
    load([vpath 'ocean_cobalt_feisty_forcing_z.',yrs,'.thkcello.mat'],'thkcello')
    thkcello = thkcello(:,:,:,1:12);

    %% Interpolate monthly forcing to daily
    ESM = daily_interp_int_monthly_means(NWID,Time,Tdays,...
        TEMP_z,TEMP_btm,det_btm,MZ,LZ,MZloss,LZloss,thkcello,WID,ni,nj);

    %% Save
    save([vpath 'ocean_cobalt_feisty_forcing_2dint.1990.mat'],'ESM','-v7.3')
end