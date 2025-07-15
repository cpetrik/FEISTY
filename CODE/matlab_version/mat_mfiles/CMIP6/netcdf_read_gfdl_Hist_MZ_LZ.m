% Read GFDL netcdfs
% CMIP6 DECK sims - Historical
% mdz & lgz biomass int over top 100m

clear
close all

fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/GFDL/hist/';

%%
ncdisp([fpath 'nmdz_100_1880_2014_all.nc'])

%%
ncdisp([fpath 'nlgz_100_1880_2014_all.nc'])

%%
ncdisp([fpath 'jhploss_nmdz_100_1880_2014_all.nc'])

%%
ncdisp([fpath 'jhploss_nlgz_100_1880_2014_all.nc'])

%%
nlgz_units      = 'molN m-2';
nmdz_units      = 'molN m-2';
lg_hp_units     = 'molN m-2 s-1';
md_hp_units     = 'molN m-2 s-1';
lz_long_name        = 'large zooplankton nitrogen biomass integral in upper 100m';
mz_long_name        = 'medium zooplankton nitrogen biomass integral in upper 100m';
lz_hploss_long_name = 'large zooplankton nitrogen loss to higher preds integral in upper 100m';
mz_hploss_long_name = 'medium zooplankton nitrogen loss to higher preds integral in upper 100m';
time_units          = 'days since 1850-01-01 00:00:00';
missing_value = 1.000000020040877e+20;

%% MZ bio
ncid = netcdf.open([fpath 'nmdz_100_1880_2014_all.nc'],...
    'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% LZ bio
ncid = netcdf.open([fpath 'nlgz_100_1880_2014_all.nc'],...
    'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% MZ loss
ncid = netcdf.open([fpath 'jhploss_nmdz_100_1880_2014_all.nc'],...
    'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% LZ loss
ncid = netcdf.open([fpath 'jhploss_nlgz_100_1880_2014_all.nc'],...
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
jhploss_nmdz_100(jhploss_nmdz_100 >= 1.00e+19) = NaN;
jhploss_nlgz_100(jhploss_nlgz_100 >= 1.00e+19) = NaN;

%%
% molN/m2 --> gC/m2
% molN/m2/s --> gC/m2/d

test1 = double(squeeze(nmdz_100(:,:,1950)))* (106/16) * 12.01;
test2 = double(squeeze(nlgz_100(:,:,1950)))* (106/16) * 12.01;
test3 = double(squeeze(jhploss_nmdz_100(:,:,1950)))* (106/16) * 12.01 *60*60*24;
test4 = double(squeeze(jhploss_nlgz_100(:,:,1950)))* (106/16) * 12.01 *60*60*24;

figure
pcolor(log10(test1)); shading flat; colorbar;
clim([0 2])

figure
pcolor(log10(test2)); shading flat; colorbar;
clim([0 2])

figure
pcolor(log10(test3)); shading flat; colorbar;
clim([-2 1])

figure
pcolor(log10(test4)); shading flat; colorbar;
clim([-2 1])

%% subset time 1950-2014
year = (time/365) + 1850;

yr65 = find(year>=1950);

mz_100 = double(nmdz_100(:,:,yr65));
lz_100 = double(nlgz_100(:,:,yr65));
hp_mz_100 = double(jhploss_nmdz_100(:,:,yr65));
hp_lz_100 = double(jhploss_nlgz_100(:,:,yr65));

%% save
save([fpath 'gfdl_hist_int100_2zmeso_monthly_1950_2014.mat'],...
    'mz_100','lz_100','hp_mz_100','hp_lz_100',...
    'mz_long_name','lz_long_name','nmdz_units','nlgz_units',...
    'mz_hploss_long_name','lz_hploss_long_name','md_hp_units','lg_hp_units',...
    'year','time','time_units','-v7.3')

% save([spath 'gfdl_hist_int100_2zmeso_monthly_1880_2014.mat'],...
%     'mdz_100','lgz_100','mz100_long_name','lz100_long_name',...
%     'nmdz_units','nlgz_units','lz100_units','mz100_units','-v7.3')

%% load and fix orientation
load([fpath 'gfdl_hist_int100_2zmeso_monthly_1950_2014.mat'])

%%
[ni,nj,nt] = size(mz_100);

% Loop over mos
mdz_100 = nan*ones(ni,nj,nt);
lgz_100 = nan*ones(ni,nj,nt);
hp_mdz_100 = nan*ones(ni,nj,nt);
hp_lgz_100 = nan*ones(ni,nj,nt);

M=0;

for t=1:nt
    testM = (squeeze(mz_100(:,:,t)));
    testL = (squeeze(lz_100(:,:,t)));

    Mtsplit1 = testM(1:180,:);
    Mtsplit2 = testM(181:end,:);
    Mtcomb = [Mtsplit2;Mtsplit1];

    Ltsplit1 = testL(1:180,:);
    Ltsplit2 = testL(181:end,:);
    Ltcomb = [Ltsplit2;Ltsplit1];

    mdz_100(:,:,t) = Mtcomb;
    lgz_100(:,:,t) = Ltcomb;

    clear Mtsplit1 Mtsplit2 Mtcomb
    clear Ltsplit1 Ltsplit2 Ltcomb


    testMH = (squeeze(hp_mz_100(:,:,t)));
    testLH = (squeeze(hp_lz_100(:,:,t)));

    Mtsplit1 = testMH(1:180,:);
    Mtsplit2 = testMH(181:end,:);
    MHcomb = [Mtsplit2;Mtsplit1];

    Ltsplit1 = testLH(1:180,:);
    Ltsplit2 = testLH(181:end,:);
    LHcomb = [Ltsplit2;Ltsplit1];

    hp_mdz_100(:,:,t) = MHcomb;
    hp_lgz_100(:,:,t) = LHcomb;

    clear Mtsplit1 Mtsplit2 MHcomb
    clear Ltsplit1 Ltsplit2 LHcomb

end

%% save
save([fpath 'gfdl_hist_int100_2zmeso_reorient_monthly_1950_2014.mat'],...
    'mdz_100','lgz_100','hp_mdz_100','hp_lgz_100',...
    'mz_long_name','lz_long_name','nmdz_units','nlgz_units',...
    'mz_hploss_long_name','lz_hploss_long_name','md_hp_units','lg_hp_units',...
    'year','time','time_units','-v7.3')



