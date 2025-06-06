% Read GFDL netcdfs
% CMIP6 DECK sims - Historical 
% mdz & lgz biomass int over top 100m

clear
close all

fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/GFDL/hist/';

%%
ncdisp([fpath 'ocean_cobalt_tracers_int.200001-200412.nmdz_100.nc'])

%%
ncdisp([fpath 'ocean_cobalt_tracers_int.200001-200412.nlgz_100.nc'])

%%
nlgz_units         = 'molN m-2';
nmdz_units         = 'molN m-2';
lz_long_name    = 'large zooplankton nitrogen biomass in upper 100m';
mz_long_name    = 'medium zooplankton nitrogen biomass in upper 100m';
lz100_units        = 'gC/m2';
mz100_units        = 'gC/m2';
lz100_long_name    = 'large zooplankton carbon biomass in upper 100m';
mz100_long_name    = 'medium zooplankton carbon biomass in upper 100m';
missing_value = 1.000000020040877e+20;

ni = 720;
nj = 576;
nt = 60;

%% loop over files, each 5 years (5*12 mo)
syr = 1850:5:2014;
eyr = 1854:5:2014;
yr = 1850:2014;
ttot = 12*length(yr); %total number of months across all files
M = 0;

mdz_100 = nan*ones(ni,nj,ttot);
lgz_100 = nan*ones(ni,nj,ttot);

for y = 1:length(yr)
    
    %%
    ncid = netcdf.open([fpath 'ocean_cobalt_tracers_int.',num2str(syr(y)),'01-',num2str(eyr(y)),'12.nmdz_100.nc'],...
        'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

    for i = 1:(nvars)
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
    end
    netcdf.close(ncid);

    %%
    ncid = netcdf.open([fpath 'ocean_cobalt_tracers_int.',num2str(syr(y)),'01-',num2str(eyr(y)),'12.nlgz_100.nc'],...
        'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

    for i = 1:(nvars)
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
    end
    netcdf.close(ncid);

    %%
    nmdz_100(nmdz_100 >= 1.00e+19) = NaN;
    nlgz_100(nlgz_100 >= 1.00e+19) = NaN;

    %% molN/kg --> gC/m3
    nmdz_100 = double(nmdz_100 * (1/1e-3) * (106/16) * 12.01);
    nlgz_100 = double(nlgz_100 * (1/1e-3) * (106/16) * 12.01);

    %% Loop over mos
    for t=1:nt
        M = M+1;

        testM = double(squeeze(nmdz_100(:,:,t)));
        testL = double(squeeze(nlgz_100(:,:,t)));

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

        mdz_100(:,:,M) = Mtcomb;
        lgz_100(:,:,M) = Ltcomb;

        clear Mtsplit1 Mtsplit2 Mtflip1 Mtflip2 Mtcomb
        clear Ltsplit1 Ltsplit2 Ltflip1 Ltflip2 Ltcomb

    end

    clear testM testL nmdz nlgz

end %yrs

%%

figure
pcolor(Mtcomb); shading flat; colorbar;
%clim([0 0.15])

figure
pcolor(Ltcomb); shading flat; colorbar;
%clim([0 0.15])

%%
test1 = double(squeeze(mdz_100(:,:,6)))./ 12.01;
test2 = double(squeeze(lgz_100(:,:,12)))./ 12.01;
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
    'mdz_100','lgz_100','mz100_long_name','lz100_long_name',...
    'nmdz_units','nlgz_units','lz100_units','mz100_units','-v7.3')





