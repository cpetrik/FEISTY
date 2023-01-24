% check netcdf save

clear 
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];

%% Setup netcdf path to store to
fname1 = 'gfdl-mom6_cobalt2_none_obsclim_histsoc_onedeg_';
fname2 = '_global_monthly_1961_2010.nc';

file_tpb = [fpath fname1 'tpb' fname2];
file_tdb = [fpath fname1 'tdb' fname2];
file_tcb = [fpath fname1 'tcb' fname2];
file_bp30 = [fpath fname1 'bp30cm' fname2];
file_bp90 = [fpath fname1 'bp90cm' fname2];
file_bd90 = [fpath fname1 'bd90cm' fname2];

file_tpc = [fpath fname1 'tpc' fname2];
file_tdc = [fpath fname1 'tdc' fname2];
file_tcc = [fpath fname1 'tcc' fname2];
file_cp30 = [fpath fname1 'cp30cm' fname2];
file_cp90 = [fpath fname1 'cp90cm' fname2];
file_cd90 = [fpath fname1 'cd90cm' fname2];

%% fb
ncdisp(file_bp30)

%% pb
ncdisp(file_tpb)

%% db
ncdisp(file_tdb)

%% tcb
ncdisp(file_tcb)

%% fc
ncdisp(file_cp30)

%% pc
ncdisp(file_tpc)

%% dc
ncdisp(file_tdc)

%% tc
ncdisp(file_tcc)

%% tcb
ncid = netcdf.open(file_tcb,'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end
netcdf.close(ncid);

%%
tcb1 = tcb(:,:,1);
tcb1(tcb1>1e19) = nan;
figure(1)
pcolor(log10(tcb1)); shading flat; colorbar
caxis([-2 3])

%% tc
ncid = netcdf.open(file_tcc,'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end
netcdf.close(ncid);

%%
tc1 = tc(:,:,1);
tc1(tc1>1e19) = nan;
figure(2)
pcolor(log10(tc1)); shading flat; colorbar
caxis([-2 3])
