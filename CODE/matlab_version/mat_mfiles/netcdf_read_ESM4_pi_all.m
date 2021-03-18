% Read GFDL ESM4 Preindust netcdfs
% all vars for FEISTY

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/ESM4_PI/';
tpath='/Volumes/MIP/GCM_DATA/ESM4_PI/Temp3D/';

%%
ncdisp([fpath 'ocean_cobalt_btm.100101-102012.fndet_btm.nc'])

% Dimensions:
% xu_ocean = 360
% yu_ocean = 200
% time     = 240   (UNLIMITED)
% nv       = 2
% xt_ocean = 360
% yt_ocean = 200
% time:units = 'days since 0001-01-01 00:00:00'

% fndet_btm
% Size:       360x200x240
% Dimensions: xt_ocean,yt_ocean,time
% Datatype:   single
% Attributes:
det_long_name     = 'ndet sinking flux to bottom';
det_units         = 'mol m-2 s-1';
% missing_value = -10000000000
% _FillValue    = -10000000000
% cell_methods  = 'time: mean'
% time_avg_info = 'average_T1,average_T2,average_DT'
% coordinates   = 'geolon_t geolat_t'

%%
ncdisp([fpath 'ocean_cobalt_btm.100101-102012.btm_temp.nc'])
% btm_temp
tb_long_name     = 'Bottom Temperature';
tb_units         = 'deg C';
%missing_value = -10000000000

%%
ncdisp([fpath 'ocean_cobalt_biomass_100.100101-102012.nmdz_100.nc'])
% nmdz_100
mz_long_name     = 'Medium zooplankton nitrogen biomass in upper 100m';
mz_units         = 'mol m-2';
%missing_value = -10000000000

%%
ncdisp([fpath 'ocean_cobalt_biomass_100.100101-102012.nlgz_100.nc'])
% nlgz_100
lz_long_name     = 'Large zooplankton nitrogen biomass in upper 100m';
lz_units         = 'mol m-2';
%missing_value = -10000000000

%%
ncdisp([fpath 'ocean_cobalt_miscflux_100.100101-102012.jhploss_nmdz_100.nc'])
% jhploss_nmdz_100
mzhp_long_name     = 'Medium zooplankton nitrogen loss to higher preds. integral in upper 100m';
mzhp_units         = 'mol m-2 s-1';
%missing_value = -10000000000

%%
ncdisp([fpath 'ocean_cobalt_miscflux_100.100101-102012.jhploss_nlgz_100.nc'])
% jhploss_nlgz_100
lzhp_long_name     = 'Large zooplankton nitrogen loss to higher preds. integral in upper 100m';
lzhp_units         = 'mol m-2 s-1';
%missing_value = -10000000000


%% Dates
%s20 = 1001:20:1481;
%e20 = 1020:20:1500; only 400 yrs of 3d temp

s20 = 1001:20:1381;
e20 = 1020:20:1400;

%% Cycle through each 20 yrs
for q = 1:length(s20)
    
    tspan = [num2str(s20(q)),'01-',num2str(e20(q)),'12'];
    
    %%
    load([tpath 'ocean.',tspan,'.temp_100.mat']);
    
    ncid = netcdf.open([fpath 'ocean_cobalt_btm.',tspan,'.fndet_btm.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
    end
    netcdf.close(ncid);
    
    ncid = netcdf.open([fpath 'ocean_cobalt_btm.',tspan,'.btm_temp.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
    end
    netcdf.close(ncid);
    
    ncid = netcdf.open([fpath 'ocean_cobalt_biomass_100.',tspan,'.nmdz_100.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
    end
    netcdf.close(ncid);
    
    ncid = netcdf.open([fpath 'ocean_cobalt_biomass_100.',tspan,'.nlgz_100.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
    end
    netcdf.close(ncid);
    
    ncid = netcdf.open([fpath 'ocean_cobalt_miscflux_100.',tspan,'.jhploss_nmdz_100.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
    end
    netcdf.close(ncid);
    
    ncid = netcdf.open([fpath 'ocean_cobalt_miscflux_100.',tspan,'.jhploss_nlgz_100.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
    end
    netcdf.close(ncid);
    
    %% NaNs on land cells
    fndet_btm(fndet_btm <= -1e9) = NaN;
    btm_temp(btm_temp <= -1e9) = NaN;
    nmdz_100(nmdz_100 <= -1e9) = NaN;
    nlgz_100(nlgz_100 <= -1e9) = NaN;
    jhploss_nmdz_100(jhploss_nmdz_100 <= -1e9) = NaN;
    jhploss_nlgz_100(jhploss_nlgz_100 <= -1e9) = NaN;
    
    %% Time
    yr = (time/365)+1; %days since 0001-01-01 = first month is Jan 1888
    t1 = round(yr(1));
    t2 = floor(yr(end));
    
    mos = length(time);
    mstart = 1:12:mos;
    mend = 12:12:mos;
    nyrs = mos/12;
    
    yrs = t1:t2;
    Tdays = 1:365;
    Time = Tdays(15:30:end);
    
    %% index of water cells
    det = fndet_btm(:,:,1);
    
    WID = find(~isnan(det(:)));
    NID = length(WID); %48111
    
    [ni,nj] = size(det);
    
    %% Make daily data
    
    % Convert units after interpolating
    %det btm: mol N m-2 s-1
    %zoo: integrated to mol N m-2
    %zoo loss: mol N m-2 s-1
    %tp: K
    %tb: degC
    
    for y = 1:nyrs
        yr = yrs(y)
        
        Tp = double(tp_100(:,:,mstart(y):mend(y)));
        Tb = double(btm_temp(:,:,mstart(y):mend(y)));
        Zm = double(nmdz_100(:,:,mstart(y):mend(y)));
        Zl = double(nlgz_100(:,:,mstart(y):mend(y)));
        dZm = double(jhploss_nmdz_100(:,:,mstart(y):mend(y)));
        dZl = double(jhploss_nlgz_100(:,:,mstart(y):mend(y)));
        det= double(fndet_btm(:,:,mstart(y):mend(y)));
        
        % setup FEISTY data files
        D_Tp  = nan*zeros(NID,365);
        D_Tb  = nan*zeros(NID,365);
        D_Zm  = nan*zeros(NID,365);
        D_Zl  = nan*zeros(NID,365);
        D_dZm  = nan*zeros(NID,365);
        D_dZl  = nan*zeros(NID,365);
        D_det = nan*zeros(NID,365);
        
        %% interpolate to daily resolution
        for j = 1:NID
            % indexes
            [m,n] = ind2sub([ni,nj],WID(j)); % spatial index of water cell
            
            % pelagic temperature (from Kelvin to Celcius)
            Y = squeeze(Tp(m,n,:)) - 273;
            yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
            D_Tp(j,:) = yi;
            
            % bottom temperature (in Celcius)
            Y = squeeze(Tb(m,n,:));
            yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
            D_Tb(j,:) = yi;
            
            % medium zoo: from mol N m-2 to g(WW) m-2
            % 106/16 mol C in 1 mol N
            % 12.01 g C in 1 mol C
            % 1 g dry W in 9 g wet W (Pauly & Christiansen)
            Y = squeeze(Zm(m,n,:));
            yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
            D_Zm(j,:) = yi * (106.0/16.0) * 12.01 * 9.0;
            
            % large zoo: from mol N m-2 to g(WW) m-2
            % 106/16 mol C in 1 mol N
            % 12.01 g C in 1 mol C
            % 1 g dry W in 9 g wet W (Pauly & Christiansen)
            Y = squeeze(Zl(m,n,:));
            yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
            D_Zl(j,:) = yi * (106.0/16.0) * 12.01 * 9.0;
            
            % medium zoo mortality: from mol N m-2 s-1 to g(WW) m-2 d-1
            % 106/16 mol C in 1 mol N
            % 12.01 g C in 1 mol C
            % 1 g dry W in 9 g wet W (Pauly & Christiansen)
            Y = squeeze(dZm(m,n,:));
            yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
            D_dZm(j,:) = yi * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
            
            % large zoo mortality: from mol N m-2 s-1 to g(WW) m-2 d-1
            % 106/16 mol C in 1 mol N
            % 12.01 g C in 1 mol C
            % 1 g dry W in 9 g wet W (Pauly & Christiansen)
            Y = squeeze(dZl(m,n,:));
            yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
            D_dZl(j,:) = yi * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
            
            % detrital flux to benthos: from mol C m-2 s-1 to g(WW) m-2 d-1
            % 106/16 mol C in 1 mol N
            % 12.01 g C in 1 mol C
            % 1 g dry W in 9 g wet W (Pauly & Christiansen)
            Y = squeeze(det(m,n,:));
            yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
            D_det(j,:) = yi * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
            
        end
        
        % Negative biomass or mortality loss from interp
        D_Zm(D_Zm<0) = 0.0;
        D_Zl(D_Zl<0) = 0.0;
        D_dZm(D_dZm<0) = 0.0;
        D_dZl(D_dZl<0) = 0.0;
        D_det(D_det<0) = 0.0;
        
        COBALT.Tp = D_Tp;
        COBALT.Tb = D_Tb;
        COBALT.Zm = D_Zm;
        COBALT.Zl = D_Zl;
        COBALT.dZm = D_dZm;
        COBALT.dZl = D_dZl;
        COBALT.det = D_det;
        
        % save
        save([fpath 'Data_ocean_cobalt_daily_',num2str(yr),'.mat'], 'COBALT');
        
    end
    
end






