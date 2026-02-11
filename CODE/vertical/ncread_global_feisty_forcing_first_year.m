% Read in first year (1st 12 mo) of forcing for feisty

function [thetao,tob,nmdz,nlgz,jhploss_n_Mdz,jhploss_n_Lgz,fntot_btm] = ncread_global_feisty_forcing_first_year(fpath,ni,nj,nz)

%%
%ncdisp([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412.nmdz.nc'])

%% thetao
ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412.thetao.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 6
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,0,0,0],[ni nj nz 12]);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

thetao(thetao>1e19) = nan;

%% tbtm
ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_2d.199001-199412.tob.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 7
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,0,0],[ni nj 12]);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

tob(tob>1e19) = nan;

%% MZ
ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412.nmdz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 5
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,0,0,0],[ni nj nz 12]);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

nmdz(nmdz>1e19) = nan;

%% LZ
ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412.nlgz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 5
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,0,0,0],[ni nj nz 12]);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

nlgz(nlgz>1e19) = nan;

%% MZ loss
ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412.jhploss_n_Mdz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 5
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,0,0,0],[ni nj nz 12]);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

jhploss_n_Mdz(jhploss_n_Mdz>1e19) = nan;

%% LZ loss
ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412.jhploss_n_Lgz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 5
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,0,0,0],[ni nj nz 12]);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

jhploss_n_Lgz(jhploss_n_Lgz>1e19) = nan;

%% Det btm
ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_2d.199001-199412.fntot_btm.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 4
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,0,0],[ni nj 12]);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

fntot_btm(fntot_btm>1e19) = nan;

end








