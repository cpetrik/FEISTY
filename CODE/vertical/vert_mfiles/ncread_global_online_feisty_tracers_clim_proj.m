% Climatol (mean of each month in 5-yr)
% MOM6-COBALTv3-FEISTY online global

clear
close all

fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';
%fpath = '/project/Feisty/NC/Global_COBALT_FEISTY/';

gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
%gpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%%
%load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet','z_l_units','z_l_long_name','z_l')

%%
ncdisp([fpath 'ocean_feisty_tracers_z.1990-1994.01.nc'])

%%
ncid = netcdf.open([fpath 'ocean_feisty_tracers_z.1990-1994.01.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 5:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'ocean_feisty_tracers_z.1990-1994.01.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 3
    varname = netcdf.inqVar(ncid, i-1);
    z_l = netcdf.getVar(ncid,i-1);
end

for i = 4
    varname = netcdf.inqVar(ncid, i-1);
    z_i = netcdf.getVar(ncid,i-1);
end
netcdf.close(ncid);

dz = diff(z_i);

%% Vertical means
vSf = mean(Sf_B,1,'omitnan');
vSf = squeeze(mean(vSf,2,'omitnan'));
vSf = squeeze(mean(vSf,2,'omitnan'));

vSp = mean(Sp_B,1,'omitnan');
vSp = squeeze(mean(vSp,2,'omitnan'));
vSp = squeeze(mean(vSp,2,'omitnan'));

vSd = mean(Sd_B,1,'omitnan');
vSd = squeeze(mean(vSd,2,'omitnan'));
vSd = squeeze(mean(vSd,2,'omitnan'));

vMf = mean(Mf_B,1,'omitnan');
vMf = squeeze(mean(vMf,2,'omitnan'));
vMf = squeeze(mean(vMf,2,'omitnan'));

vMp = mean(Mp_B,1,'omitnan');
vMp = squeeze(mean(vMp,2,'omitnan'));
vMp = squeeze(mean(vMp,2,'omitnan'));

vLp = mean(Lp_B,1,'omitnan');
vLp = squeeze(mean(vLp,2,'omitnan'));
vLp = squeeze(mean(vLp,2,'omitnan'));

vLd = mean(Ld_B,1,'omitnan');
vLd = squeeze(mean(vLd,2,'omitnan'));
vLd = squeeze(sum(vLd,2,'omitnan'));


%% vert sums
iSf = squeeze(sum((Sf_B.*thkcello),3,'omitnan'));
iSp = squeeze(sum((Sp_B.*thkcello),3,'omitnan'));
iSd = squeeze(sum((Sd_B.*thkcello),3,'omitnan'));
iMf = squeeze(sum((Mf_B.*thkcello),3,'omitnan'));
iMp = squeeze(sum((Mp_B.*thkcello),3,'omitnan'));
iLp = squeeze(sum((Lp_B.*thkcello),3,'omitnan'));

%% spatial mean (of vert integral)
sSf = mean(iSf,3,'omitnan');
sSp = mean(iSp,3,'omitnan');
sSd = mean(iSd,3,'omitnan');
sMf = mean(iMf,3,'omitnan');
sMp = mean(iMp,3,'omitnan');
sLp = mean(iLp,3,'omitnan');

sMd = mean(Md_B(:,:,35,:),4,'omitnan');
sLd = mean(Ld_B(:,:,35,:),4,'omitnan');
sBe = mean(BE_B(:,:,35,:),4,'omitnan');

