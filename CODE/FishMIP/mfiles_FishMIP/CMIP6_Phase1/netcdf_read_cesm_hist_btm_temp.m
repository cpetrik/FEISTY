% Read Fish-MIP CESM Historic netcdfs
% Concatenate bottom temperature

clear all
close all

fpath='/Volumes/GFDL/Fish-MIP/CESM/Hist/';

tstart = 185001:1000:200512;
tend = 185912:1000:200512;
tend = [tend 200512];

%each file has 10 years of 12 months, except last is 6 yrs of 12 months
mos = 10*12*(length(tstart) - 1) + 6*12;
mod_time = nan(mos,1);
tbnds = nan(2,mos);
nsmz_100 = nan(360,180,mos);

mstart = 1:120:mos;
mend = [120:120:mos mos];

%%
for t=1:length(tstart)
ncid = netcdf.open([fpath 'cesm_hist_to_zb_monthly_',num2str(tstart(t)),'-',...
    num2str(tend(t)),'.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%NaNs on land cells
to(to>=9e+36)=nan;
tb = to;
save([fpath 'cesm_hist_tb_100_monthly_',num2str(tstart(t)),'-',...
    num2str(tend(t)),'.mat'],'tb','time');

%concatenate
mod_time(mstart(t):mend(t)) = time;
% tbnds(:,mstart(t):mend(t)) = time_bnds;
tbtm(:,:,mstart(t):mend(t)) = tb;

clear to time time_bnds 

end

save([fpath 'cesm_hist_tbtm_monthly_',num2str(tstart(1)),'-',...
    num2str(tend(end)),'.mat'],'tbtm','mod_time');
