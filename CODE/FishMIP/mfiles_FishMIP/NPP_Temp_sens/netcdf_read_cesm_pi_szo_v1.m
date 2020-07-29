% Read Fish-MIP netcdfs
% Integrate over top 100 m for small zooplankton

clear all
close all

%%
fpath='/Volumes/GFDL/Fish-MIP/CESM/PreIndust/';
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_185001-185912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

szoo(szoo>=9e+36)=nan;
time1 = time;
tbnds1 = time_bnds;
%top 100 m
sz100 = szoo(:,:,1:10,:);
%integrated over top 100m
isz100 = squeeze(nansum(sz100,3));

save([fpath 'cesm_pi_szoo_100_185001-185912.m'],'isz100','time','time_bnds');

clear szoo time time_bnds 

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_186001-186912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_187001-187912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_188001-188912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_189001-189912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_190001-190912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_191001-191912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_192001-192912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_193001-193912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_194001-194912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_195001-195912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_196001-196912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_197001-197912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_198001-198912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_199001-199912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_200001-200912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_201001-201912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_202001-202912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_203001-203912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_204001-204912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_205001-205912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_206001-206912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_207001-207912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_208001-208912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_209001-209912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
ncid = netcdf.open([fpath 'cesm_pi_szoo_zall_monthly_210001-210012.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);


%% Calculate quantities
dat = dxt.*dyt;
datr = 1.0./(dat+eps);

tmask = wet;

%% Save needed variables
save([fpath fname '.mat'])

%% Retain all 2D
GRD.LON = double(geolon);
GRD.LAT = double(geolat);
GRD.Z   = double(deptho);
GRD.DX = double(dxt);
GRD.DY = double(dyt);
GRD.AREA  = double(areacello);
GRD.dxtn  = double(dxt);
GRD.dyte  = double(dyt);
GRD.datr  = double(datr);
GRD.lmask = double(tmask);

%! save
save([fpath 'Data_grid2D_MOM6_preindust.mat'],'GRD');

%% Retain 2D for mom6 diffusion
[ni,nj] = size(dat);
G.isc = 1;
G.iec = ni;
G.jsc = 2; %1st grid cell is Antarctica
G.jec = nj;
%don't know what these are, assume same as above
G.IscB = 1;
G.IecB = ni;
G.JscB = 2;
G.JecB = nj;
G.dxCu = double(dxCu);
G.dxCv = double(dxCv);
G.dyCu = double(dyCu);
G.dyCv = double(dyCv);
G.area = double(areacello);
G.Iarea = double(datr);
G.dxt  = double(dxt);
G.dyt  = double(dyt);
G.mask = double(tmask);

%! save
save([fpath 'Data_grid2D_diff_MOM6_preindust.mat'],'G');

%% Retain only water cells and vectorize
clear GRD
ID = find(wet(:)>0);
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = geolon(ID);
GRD.LAT = geolat(ID);
GRD.Z   = deptho(ID);
GRD.DX = dxt(ID);
GRD.DY = dyt(ID);
GRD.AREA  = areacello(ID);
GRD.dxtn  = dxt(ID);
GRD.dyte  = dyt(ID);
GRD.datr  = datr(ID);
GRD.lmask = tmask(ID);

%! save
save([fpath 'Data_grid1D_MOM6_preindust.mat'],'GRD');
        



