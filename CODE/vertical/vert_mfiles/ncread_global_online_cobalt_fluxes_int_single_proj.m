% 0.5 degree global MOM6-COBALTv3-FEISTY only
% C, N, O vars

clear
close all

%fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';
%fpath = '/project/Feisty/NC/Global_COBALT_FEISTY/cobalt_feisty/';
fpath = '/project/Feisty/Globus_RW/COBALT-FEISTY/';

%gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
gpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

spath = '/project/Feisty/NC/Global_COBALT_FEISTY/cobalt_feisty/';

%%
%load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet','z_l_units','z_l_long_name','z_l')

%%
ncdisp([fpath '19900101.ocean_cobalt_fluxes_int.nc'])

%%
ncid = netcdf.open([fpath '19900101.ocean_cobalt_fluxes_int.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% 
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

% wc_vert_int_poc(wc_vert_int_poc>1e19) = nan;
% wc_vert_int_doc(wc_vert_int_doc>1e19) = nan;
% wc_vert_int_dic(wc_vert_int_dic>1e19) = nan;
% wc_vert_int_o2(wc_vert_int_o2>1e19) = nan;

wc_vert_int_npp(wc_vert_int_npp>1e19) = nan;
jzloss_nmdz_100(jzloss_nmdz_100>1e19) = nan;
jprod_nmdz_100(jprod_nmdz_100>1e19) = nan;
jprod_nlgz_100(jprod_nlgz_100>1e19) = nan;

jingest_n_hp_100(jingest_n_hp_100>1e19) = nan; 
jhploss_nmdz_100(jhploss_nmdz_100>1e19) = nan;
jhploss_nlgz_100(jhploss_nlgz_100>1e19) = nan;

%% Time series of vert integral
tNPP = mean(wc_vert_int_npp,1,'omitnan');
tNPP = squeeze(mean(tNPP,2,'omitnan'));

tMZzl = mean(jzloss_nmdz_100,1,'omitnan');
tMZzl = squeeze(mean(tMZzl,2,'omitnan'));

tMZprod = mean(jprod_nmdz_100,1,'omitnan');
tMZprod = squeeze(mean(tMZprod,2,'omitnan'));

tLZprod = mean(jprod_nlgz_100,1,'omitnan');
tLZprod = squeeze(mean(tLZprod,2,'omitnan'));

tHPinj = mean(jingest_n_hp_100,1,'omitnan');
tHPinj = squeeze(mean(tHPinj,2,'omitnan'));

tMZloss = mean(jhploss_nmdz_100,1,'omitnan');
tMZloss = squeeze(mean(tMZloss,2,'omitnan'));

tLZloss = mean(jhploss_nlgz_100,1,'omitnan');
tLZloss = squeeze(mean(tLZloss,2,'omitnan'));

%% spatial mean of vert integral
sNPP = mean(wc_vert_int_npp,3,'omitnan');
sMZzl = mean(jzloss_nmdz_100,3,'omitnan');
sMZprod = mean(jprod_nmdz_100,3,'omitnan');
sLZprod = mean(jprod_nlgz_100,3,'omitnan');

sHPinj = mean(jingest_n_hp_100,3,'omitnan');
sMZloss = mean(jhploss_nmdz_100,3,'omitnan');
sLZloss = mean(jhploss_nlgz_100,3,'omitnan');

%%
save([spath '19900101.ocean_cobalt_fluxes_int_means.mat'],...
    'tNPP','tMZzl','tMZprod','tLZprod',...
    'sNPP','sMZzl','sMZprod','sLZprod',...
    'tHPinj','tMZloss','tLZloss',...
    'sHPinj','sMZloss','sLZloss')






