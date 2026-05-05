% 0.5 degree global MOM6-COBALTv3-FEISTY online
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
ncdisp([fpath '19940101.ocean_cobalt_tracers_month_z.nc'])

%%
%
ncid = netcdf.open([fpath '19940101.ocean_cobalt_tracers_month_z.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% xh and yh
for i = 1:2
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end

% time
for i = 5
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end

% phyto and zoo
for i = 8:14
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end

% ndet
for i = 25
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end

% c,p,o
for i = 39:43
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

ndet(ndet>1e19) = nan;
nsm(nsm>1e19) = nan;
nlg(nlg>1e19) = nan;
ndi(ndi>1e19) = nan;
nsmz(nsmz>1e19) = nan;
no3(no3>1e19) = nan;
nh4(nh4>1e19) = nan;
chl(chl>1e19) = nan;
schl = squeeze(chl(:,:,1,:));
o2(o2>1e19) = nan;

%% thkcello
load([gpath 'ocean_cobalt_feisty_forcing_z.199001-199412.thkcello.mat'])
%thkcello = thkcello(:,:,:,1:12);
thkcello = thkcello(:,:,:,(60-11):60);

%load([gpath 'ocean_cobalt_feisty_forcing_z.199501-199912.thkcello.mat'])


%% Vertical means
vNo3 = mean(no3,1,'omitnan');
vNo3 = squeeze(mean(vNo3,2,'omitnan'));
vNo3 = squeeze(mean(vNo3,2,'omitnan'));

vNh4 = mean(nh4,1,'omitnan');
vNh4 = squeeze(mean(vNh4,2,'omitnan'));
vNh4 = squeeze(mean(vNh4,2,'omitnan'));

vChl = mean(chl,1,'omitnan');
vChl = squeeze(mean(vChl,2,'omitnan'));
vChl = squeeze(mean(vChl,2,'omitnan'));

vo2 = mean(o2,1,'omitnan');
vo2 = squeeze(mean(vo2,2,'omitnan'));
vo2 = squeeze(mean(vo2,2,'omitnan'));

vSp = mean(nsm,1,'omitnan');
vSp = squeeze(mean(vSp,2,'omitnan'));
vSp = squeeze(mean(vSp,2,'omitnan'));

vLp = mean(nlg,1,'omitnan');
vLp = squeeze(mean(vLp,2,'omitnan'));
vLp = squeeze(mean(vLp,2,'omitnan'));

vSz = mean(nsmz,1,'omitnan');
vSz = squeeze(mean(vSz,2,'omitnan'));
vSz = squeeze(mean(vSz,2,'omitnan'));

vDe = mean(ndet,1,'omitnan');
vDe = squeeze(mean(vDe,2,'omitnan'));
vDe = squeeze(mean(vDe,2,'omitnan'));

vDi = mean(ndi,1,'omitnan');
vDi = squeeze(mean(vDi,2,'omitnan'));
vDi = squeeze(mean(vDi,2,'omitnan'));

%% vert sums or means
mNo3 = squeeze(sum((no3.*thkcello),3,'omitnan')) ./ squeeze(sum(thkcello,3,'omitnan'));
mNh4 = squeeze(sum((nh4.*thkcello),3,'omitnan')) ./ squeeze(sum(thkcello,3,'omitnan'));
iChl = squeeze(sum((chl.*thkcello),3,'omitnan'));
mo2 = squeeze(sum((o2.*thkcello),3,'omitnan')) ./ squeeze(sum(thkcello,3,'omitnan'));

iSp = squeeze(sum((nsm.*thkcello),3,'omitnan'));
iLp = squeeze(sum((nlg.*thkcello),3,'omitnan'));
iSz = squeeze(sum((nsmz.*thkcello),3,'omitnan'));
iDe = squeeze(sum((ndet.*thkcello),3,'omitnan'));
iDi = squeeze(sum((ndi.*thkcello),3,'omitnan'));

%% Time series of vert integral
tNo3 = mean(mNo3,1,'omitnan');
tNo3 = squeeze(mean(tNo3,2,'omitnan'));

tNh4 = mean(mNh4,1,'omitnan');
tNh4 = squeeze(mean(tNh4,2,'omitnan'));

tChl = mean(iChl,1,'omitnan');
tChl = squeeze(mean(tChl,2,'omitnan'));

to2 = mean(mo2,1,'omitnan');
to2 = squeeze(mean(to2,2,'omitnan'));

tSp = mean(iSp,1,'omitnan');
tSp = squeeze(mean(tSp,2,'omitnan'));

tLp = mean(iLp,1,'omitnan');
tLp = squeeze(mean(tLp,2,'omitnan'));

tSz = mean(iSz,1,'omitnan');
tSz = squeeze(mean(tSz,2,'omitnan'));

tDe = mean(iDe,1,'omitnan');
tDe = squeeze(mean(tDe,2,'omitnan'));

tDi = mean(iDi,1,'omitnan');
tDi = squeeze(mean(tDi,2,'omitnan'));

%% spatial mean of vert integral
sNo3 = mean(mNo3,3,'omitnan');
sNh4 = mean(mNh4,3,'omitnan');
so2  = mean(mo2,3,'omitnan');
sChls = mean(schl,3,'omitnan');

sChl = mean(iChl,3,'omitnan');
sSp = mean(iSp,3,'omitnan');
sLp = mean(iLp,3,'omitnan');
sSz = mean(iSz,3,'omitnan');
sDe = mean(iDe,3,'omitnan');
sDi = mean(iDi,3,'omitnan');


%%
save([spath '19940101.ocean_cobalt_tracers_month_z.mat'],...
    'tNo3','tNh4','tChl','to2',...
    'sNo3','sNh4','sChl','sChls','so2',...
    'vNo3','vNh4','vChl','vo2',...
    'tSp','tLp','tSz','tDe','tDi',...
    'sSp','sLp','sSz','sDe','sDi',...
    'vSp','vLp','vSz','vDe','vDi')






