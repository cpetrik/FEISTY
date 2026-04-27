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
%ncdisp([fpath '19900101.ocean_cobalt_feisty_forcing_z.nc'])
%ncdisp([fpath '19900101.ocean_cobalt_feisty_forcing_2d.nc'])

%%
nt=12;

%% 3D
ncid = netcdf.open([fpath '19900101.ocean_cobalt_feisty_forcing_z.nc'],'NC_NOWRITE');
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

thetao(thetao>1e19) = nan;
nmdz(nmdz>1e19) = nan;
nlgz(nlgz>1e19) = nan;
jhploss_n_Mdz(jhploss_n_Mdz>1e19) = nan;
jhploss_n_Lgz(jhploss_n_Lgz>1e19) = nan;

%% Btm
ncid = netcdf.open([fpath '19900101.ocean_cobalt_feisty_forcing_2d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 5:6
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

tob(tob>1e19) = nan;
fntot_btm(fntot_btm>1e19) = nan;

%% thkcello
load([gpath 'ocean_cobalt_feisty_forcing_z.199001-199412.thkcello.mat'])
thkcello = thkcello(:,:,:,1:12);

%% Vertical means
vMz = mean(nmdz,1,'omitnan');
vMz = squeeze(mean(vMz,2,'omitnan'));
vMz = squeeze(mean(vMz,2,'omitnan'));

vMhp = mean(jhploss_n_Mdz,1,'omitnan');
vMhp = squeeze(mean(vMhp,2,'omitnan'));
vMhp = squeeze(mean(vMhp,2,'omitnan'));

vLz = mean(nlgz,1,'omitnan');
vLz = squeeze(mean(vLz,2,'omitnan'));
vLz = squeeze(mean(vLz,2,'omitnan'));

vLhp = mean(jhploss_n_Lgz,1,'omitnan');
vLhp = squeeze(mean(vLhp,2,'omitnan'));
vLhp = squeeze(mean(vLhp,2,'omitnan'));

%% vert sums or means
vMZ = squeeze(sum((nmdz.*thkcello),3,'omitnan'));
vMH = squeeze(sum((jhploss_n_Mdz.*thkcello),3,'omitnan'));

vLZ = squeeze(sum((nlgz.*thkcello),3,'omitnan'));
vLH = squeeze(sum((jhploss_n_Lgz.*thkcello),3,'omitnan'));

mTP = squeeze(sum((thetao(:,:,1:10,:).*thkcello(:,:,1:10,:)),3,'omitnan') ./ sum(thkcello(:,:,1:10,:),3,'omitnan'));

%% Time series of vert integral
tMz = mean(vMZ,1,'omitnan');
tMz = squeeze(mean(tMz,2,'omitnan'));

tMhp = mean(vMH,1,'omitnan');
tMhp = squeeze(mean(tMhp,2,'omitnan'));

tLz = mean(vLZ,1,'omitnan');
tLz = squeeze(mean(tLz,2,'omitnan'));

tLhp = mean(vLH,1,'omitnan');
tLhp = squeeze(mean(tLhp,2,'omitnan'));

tTp = mean(mTP,1,'omitnan');
tTp = squeeze(mean(tTp,2,'omitnan'));

tTb = mean(tob,1,'omitnan');
tTb = squeeze(mean(tTb,2,'omitnan'));

tDet = mean(fntot_btm,1,'omitnan');
tDet = squeeze(mean(tDet,2,'omitnan'));

%% spatial mean of vert integral
sMz = mean(vMZ,3,'omitnan');
sMhp = mean(vMH,3,'omitnan');
sLz = mean(vLZ,3,'omitnan');
sLhp = mean(vLH,3,'omitnan');
sTp = mean(mTP,3,'omitnan');
sTb = mean(tob,3,'omitnan');
sDet = mean(fntot_btm,3,'omitnan');

%%
save([fpath 'ocean_cobalt_feisty_forcing_1990_means.mat'],...
    'tMz','tMhp','tLz','tLhp','tTp','tTb','tDet',...
    'sMz','sMhp','sLz','sLhp','sTp','sTb','sDet',...
    'vMZ','vMH','vLZ','vLH','mTP')






