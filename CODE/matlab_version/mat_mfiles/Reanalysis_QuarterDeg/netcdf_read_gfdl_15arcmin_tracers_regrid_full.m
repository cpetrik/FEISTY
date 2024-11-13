% Read GFDL netcdfs
% obsclim
% mdz & lgz biomass all depths
% Regridded by code merging Xiao's with 1 line from Matthias

clear
close all

fpath='/Volumes/petrik-lab/jabrzenski/Remapping_netCDF/FLUX/';
qpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
spath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';

%%
nlgz_units         = 'molN/kg';
nmdz_units         = 'molN/kg';
lz100_units        = 'gC/m2';
mz100_units        = 'gC/m2';
lz100_long_name    = 'large zooplankton biomass integrated over top 100m';
mz100_long_name    = 'medium zooplankton biomass integrated over top 100m';
missing_value = 1.000000020040877e+20;

%% Cell thickness
ni = 1440;
nj = 720;
nt = 12;

load([qpath 'gfdl-mom6-cobalt2_ctrlclim_thkcello_15arcmin_global_fixed.mat'])

% thkcello is the same everywhere, can take mean
thkcello(thkcello >= 1.00e+20) = NaN;
thk = (mean(thkcello,1,'omitnan'));
thk = (mean(thk,2,'omitnan'));

thk_mat = repmat(thk,ni,nj,1,nt);
%thk_100 = thk_mat(:,:,z100,:);

%% loop over files, each one year (12 mo)
yr = 1961:2010;
ttot = 12*length(yr); %total number of months across all files
M = 0;

nmdz_100 = nan*ones(ni,nj,ttot);
nlgz_100 = nan*ones(ni,nj,ttot);

for y = 1:length(yr)
    Y = yr(y);

    %%
    ncid = netcdf.open([fpath num2str(Y) '0101.ocean_cobalt_tracers_month_z_FishMIP_remapped.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

    for i = 1:(nvars-2)
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
    end

    %% Get top 100 m
    z100 = find(z_l <= 100);
    nk = length(z_l);
    thk_100 = thk_mat(:,:,z100,:);

    %%
    for i = (nvars-1):nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1, [0,0,0,0],[ni nj length(z100) nt]);']);
        eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
    end
    netcdf.close(ncid);

    %%
    nmdz(nmdz >= 1.00e+19) = NaN;
    nlgz(nlgz >= 1.00e+19) = NaN;

    %% molN/kg --> gC/m3
    nmdz = double(nmdz * (1/1e-3) * (106/16) * 12.01);
    nlgz = double(nlgz * (1/1e-3) * (106/16) * 12.01);

    %% Integrate top 100 m
    mz100 = double(squeeze(sum(nmdz .* thk_100 , 3, 'omitnan')));
    lz100 = double(squeeze(sum(nlgz .* thk_100 , 3, 'omitnan')));

    if y == 1
        omask = squeeze(mz100(:,:,6));
        omask(omask<=0) = nan;
        omask(omask>0) = 1;
        omask = repmat(omask,1,1,nt);
    end

    mz100 = mz100 .* omask;
    lz100 = lz100 .* omask;

    %% Loop over mos
    for t=1:nt
        M = M+1;

        testM = double(squeeze(mz100(:,:,t)));
        testL = double(squeeze(lz100(:,:,t)));

        Mtsplit1 = testM(1:720,:);
        Mtsplit2 = testM(721:end,:);
        Mtflip1 = fliplr(Mtsplit1);
        Mtflip2 = fliplr(Mtsplit2);
        Mtcomb = [Mtflip2;Mtflip1];

        Ltsplit1 = testL(1:720,:);
        Ltsplit2 = testL(721:end,:);
        Ltflip1 = fliplr(Ltsplit1);
        Ltflip2 = fliplr(Ltsplit2);
        Ltcomb = [Ltflip2;Ltflip1];

        nmdz_100(:,:,M) = Mtcomb;
        nlgz_100(:,:,M) = Ltcomb;

        clear Mtsplit1 Mtsplit2 Mtflip1 Mtflip2 Mtcomb
        clear Ltsplit1 Ltsplit2 Ltflip1 Ltflip2 Ltcomb

    end

    clear mz100 lz100 nmdz nlgz

end %yrs

%%
test1 = double(squeeze(nmdz_100(:,:,6)))./ 12.01;
test2 = double(squeeze(nlgz_100(:,:,12)))./ 12.01;
test3 = double(squeeze(thkcello(:,:,5)));

figure
pcolor(test1); shading flat; colorbar;
clim([0 0.15])

figure
pcolor(test2); shading flat; colorbar;
clim([0 0.15])

figure
pcolor(test3); shading flat

%% save
save([spath '19610101-20101231.ocean_cobalt_tracers_int100_FishMIP_remapped.mat'],...
    'nmdz_100','nlgz_100','mz100_long_name','lz100_long_name',...
    'nmdz_units','nlgz_units','lz100_units','mz100_units','-v7.3')





