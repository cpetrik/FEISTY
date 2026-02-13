% 0.5 degree global MOM6-COBALTv3 only


clear
close all

fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%% years
ystart = 1990:5:2015;
yend = 1994:5:2019;
nyears = length(ystart);

%%
%ncdisp([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412.volcello.nc'])

%% area
load([fpath 'grid_OM4_05_COBALTv3.mat'],'areacello','areacello_long_name','areacello_units');

%% loop over years
for y = 5:6 %1:nyears

    %% volume
    ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_z.',num2str(ystart(y)),'01-',num2str(yend(y)),'12.volcello.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 8 %1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    volcello(volcello>1e19) = nan;

    %% calc thk

    [ni,nj,nk,nt] = size(volcello);
    area_mat = repmat(areacello,1,1,nk,nt);

    thkcello = volcello./area_mat;

    %% save

    save([fpath 'ocean_cobalt_feisty_forcing_z.',num2str(ystart(y)),'01-',num2str(yend(y)),'12.thkcello.mat'],...
        'thkcello',"-v7.3")


end
