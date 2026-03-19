% 0.5 degree global MOM6-COBALTv3-FEISTY online


clear
close all

%fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/cobalt_feisty/';
fpath = '/project/Feisty/NC/Global_COBALT_FEISTY/cobalt_feisty/';

%gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
gpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%%
%load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet','z_l_units','z_l_long_name','z_l')

%%
%ncdisp([fpath '19900101.ocean_feisty_pelagic_fluxes_z.nc'])

%% want met, prod, E_A, f_tot, Fout, rho

ncid = netcdf.open([fpath '19900101.ocean_feisty_forage_fluxes_z.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:2
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
for i = 5:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

ncid = netcdf.open([fpath '19900101.ocean_feisty_pelagic_fluxes_z.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 5:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

ncid = netcdf.open([fpath '19900101.ocean_feisty_demersal_fluxes_z.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 5:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

%% thkcello
load([gpath 'ocean_cobalt_feisty_forcing_z.199001-199412.thkcello.mat'])

%% 
% Define the fish groups and the base variables
fishTypes = {'SF', 'MF', 'SP', 'MP', 'LP', 'SD', 'MD', 'LD', 'BE'};
baseVars = {'met', 'prod', 'E_A', 'f_tot', 'Fout', 'rho'};

% Initialize structures to store the processed outputs
ts_means = struct();        % Will hold [time, 1] arrays
spat_vert_ints = struct();  % Will hold [lon, lat] arrays
vert_means = struct();      % Will hold [depth, 1] arrays

% Loop through each fish type
for f = 1:length(fishTypes)
    fish = fishTypes{f};
    
    % Loop through each base variable for that fish type
    for v = 1:length(baseVars)
        base = baseVars{v};
        
        % Construct the exact variable name (e.g., 'SF_met')
        varName = [fish, '_', base];
        
        % Fetch the 4D array from the base workspace using the string name
        % (Make sure thkcello is also loaded in your workspace)
        try
            currentVar = evalin('base', varName);
            currentVar(currentVar > 1e19) = nan;
        catch
            warning('Variable %s not found in workspace. Skipping...', varName);
            continue; % Skip to the next iteration if the variable is missing
        end
        
        %% 1. Time series means [time, 1]
        % Average over lon (1), lat (2), and depth (3). 
        tmp_ts = squeeze(mean(currentVar, [1, 2, 3], 'omitnan'));
        ts_means.(varName) = tmp_ts(:);
        
        %% 2. Spatial vertical integrals averaged over time [lon, lat]
        % Multiply by thkcello, sum over depth (3), then average over time (4)
        int_depth = sum(currentVar .* thkcello, 3, 'omitnan');
        spat_vert_ints.(varName) = squeeze(mean(int_depth, 4, 'omitnan')); 
        
        %% 3. Vertical means averaged over all time and horizontal space [depth, 1]
        % Average over lon (1), lat (2), and time (4).
        tmp_vert = squeeze(mean(currentVar, [1, 2, 4], 'omitnan'));
        vert_means.(varName) = tmp_vert(:);
        
    end
end

%%
save([fpath '19900101.ocean_feisty_fluxes_z_means.mat'],...
    'fishTypes','baseVars','ts_means','spat_vert_ints','vert_means')






