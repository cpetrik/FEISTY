% Read Fish-MIP RCP 8.5 netcdfs
% Mean top 100 m for temperature

clear all
close all

fpath='/Volumes/GFDL/Fish-MIP/CESM/RCP85/';

tstart = 201001:1000:210001;
tend = 200912:1000:210012;
tstart = [200601 tstart];
tend = [tend 210012];

%2006-2010 is 4 years, then
%each file has 10 years of 12 months, 
%last is 1 year of 12 months
mos = 4*12 + 10*12*(length(tstart) - 2) + 12;
mod_time = nan(mos,1);
tbnds = nan(2,mos);
nsmz_100 = nan(360,180,mos);

mstart = [1 (1+48):120:mos];
mend = [48:120:mos mos];

%%
for t=1:length(tstart)
ncid = netcdf.open([fpath 'cesm_rcp85_to_zall_monthly_',num2str(tstart(t)),'-',...
    num2str(tend(t)),'.nc4'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%NaNs on land cells
to(to>=9e+36)=nan;
%top 100 m
tp100 = to(:,:,1:10,:);
%integrated over top 100m
itp100 = squeeze(nanmean(tp100,3));

save([fpath 'cesm_rcp85_tp_100_monthly_',num2str(tstart(t)),'-',...
    num2str(tend(t)),'.mat'],'itp100','time','time_bnds');

%concatenate
mod_time(mstart(t):mend(t)) = time;
tbnds(:,mstart(t):mend(t)) = time_bnds;
tp_100(:,:,mstart(t):mend(t)) = itp100;

clear to time time_bnds 

end

save([fpath 'cesm_rcp85_tp_100_monthly_',num2str(tstart(1)),'-',...
    num2str(tend(end)),'.mat'],'tp_100','mod_time','tbnds');
