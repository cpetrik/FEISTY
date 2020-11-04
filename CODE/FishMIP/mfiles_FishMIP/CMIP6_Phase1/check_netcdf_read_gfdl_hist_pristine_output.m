% FEISTY output at all locations

clear all
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/MIP/NC/FishMIP/GFDL_CMIP6/' cfile '/'];

%% Setup netcdf path to store to
fname1 = 'feisty_gfdl-esm4_nobasd_historical_nat_default_';
fname2 = '_global_monthly_1950_2014.nc';

file_tpb = [fpath fname1 'tpb' fname2];
file_tdb = [fpath fname1 'tdb' fname2];
file_tcb = [fpath fname1 'tcb' fname2];
file_bp30 = [fpath fname1 'bp30cm' fname2];
% file_bp3090 = [fpath fname1 'bp30to90cm' fname2];
file_bp90 = [fpath fname1 'bp90cm' fname2];
% file_bd30 = [fpath fname1 'bd30cm' fname2];
% file_bd3090 = [fpath fname1 'bd30to90cm' fname2];
file_bd90 = [fpath fname1 'bd90cm' fname2];

%% 
ncdisp(file_tpb)
                         
%% 
ncid = netcdf.open(file_tcb,'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% map
figure
axesm ('mollweid','MapLatLimit',[-90 90],'MapLonLimit',[-270 90],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat,lon,log10(tcb(:,:,700)))
colormap('jet')
colorbar
caxis([-1 2])
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% write example csvs
writematrix(LAT,[fpath fname1 'LAT_global_monthly_1950-2014.csv'],...
    'Delimiter',',')

writematrix(LON,[fpath fname1 'LON_global_monthly_1950-2014.csv'],...
    'Delimiter',',')

tcb700 = tcb(:,:,700);
writematrix(tcb700,[fpath fname1 'tcb700_global_monthly_1950-2014.csv'],...
    'Delimiter',',')




