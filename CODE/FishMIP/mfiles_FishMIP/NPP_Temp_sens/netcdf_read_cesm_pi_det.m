% Read Fish-MIP netcdfs
% Concatenate bottom POC flux
% "POC Flux into Cell on Bottom (z_b)" 
% POC_FLUX_IN:units = "mmol/m^3 cm/s"

clear all
close all

fpath='/Volumes/GFDL/Fish-MIP/CESM/PreIndust/';

tstart = 185001:1000:210001;
tend = 185912:1000:210912;
tend(end)=210012;

%each file has 10 years of 12 months, except last is 1 year of 12 months
mos = 10*12*(length(tstart) - 1) + 12;
mod_time = nan(mos,1);
%tbnds = nan(2,mos);
poc_btm = nan(360,180,mos);

mstart = 1:120:mos;
mend = [120:120:mos 3012];

%%
for t=1:length(tstart)
ncid = netcdf.open([fpath 'cesm_pi_POC_FLUX_IN_zb_monthly_',num2str(tstart(t)),'-',...
    num2str(tend(t)),'.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%NaNs on land cells
det = POC_FLUX_IN;
save([fpath 'cesm_pi_det_btm_monthly_',num2str(tstart(t)),'-',...
    num2str(tend(t)),'.mat'],'det','time');

%concatenate
mod_time(mstart(t):mend(t)) = time;
% tbnds(:,mstart(t):mend(t)) = time_bnds;
poc_btm(:,:,mstart(t):mend(t)) = det;

clear POC_FLUX_IN time 

end

save([fpath 'cesm_pi_det_btm_monthly_',num2str(tstart(1)),'-',...
    num2str(tend(end)),'.mat'],'poc_btm','mod_time');




