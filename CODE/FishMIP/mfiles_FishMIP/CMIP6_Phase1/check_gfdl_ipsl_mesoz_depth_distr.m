% Look at depth distribution of snipit of zoop
% Across IPSL and GFDL
% Have similar zmeso vint, but GFDL more in <100m

clear all
close all

%% GFDL
fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
 
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_historical_zmeso_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Time
yr = ((time+1)/12)+1601;
runs = find(yr>2000 & yr<=2001);
z100 = lev;

%% Get subset of time & depth
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs(1)-1],[360 180 length(z100) length(runs)]);
    
end
netcdf.close(ncid);
zmeso(zmeso >= 1.00e+20) = NaN;

%%
gzoo = squeeze(nanmean(nanmean(nanmean(zmeso,4),2),1));
%%
plot(log10(gzoo),-lev)
clear zmeso

%% IPSL
ipath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
 
ncid = netcdf.open([ipath 'ipsl-cm6a-lr_r1i1p1f1_historical_zmeso_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Time
yr = ((time+1)/12)+1601;
runs = find(yr>2000 & yr<=2001);
z100 = olevel;

%% Get subset of time & depth
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs(1)-1], [360 180 length(z100) length(runs)]);
end
netcdf.close(ncid);
zmeso(zmeso >= 1.00e+20) = NaN;

%%
whos zmeso
izoo = squeeze(nanmean(nanmean(nanmean(zmeso,4),2),1));

%%
plot(log10(izoo),-olevel)
clear zmeso

%%
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP6/';

figure(3)
plot(log10(gzoo),-lev,'b'); hold on
plot(log10(izoo),-olevel,'r')
ylim([-500 0])
ylabel('Depth')
xlabel('log_1_0 mol C m^-^3')
title('Zmeso in Hist 2000')
legend('GFDL','IPSL')
legend('location','northwest')
print('-dpng',[pp 'Zmeso_log_depth_distr_2000.png'])

%%
figure(4)
plot((gzoo),-lev,'b'); hold on
plot((izoo),-olevel,'r')
ylim([-500 0])
ylabel('Depth')
xlabel('mol C m^-^3')
title('Zmeso in Hist 2000')
legend('GFDL','IPSL')
legend('location','northwest')
print('-dpng',[pp 'Zmeso_depth_distr_2000.png'])

%%



