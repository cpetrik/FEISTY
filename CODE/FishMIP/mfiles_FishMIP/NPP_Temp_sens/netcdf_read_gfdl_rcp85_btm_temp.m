% Read Fish-MIP RCP 8.5 netcdfs
% Concatenate bottom temp

clear all
close all

fpath='/Volumes/FEISTY/Fish-MIP/GFDL/RCP85/';
cpath='/Volumes/GFDL/Fish-MIP/GFDL/RCP85/';

tstart = 201101:1000:210000;
tend = 201012:1000:210012;
tstart = [200601 tstart];

%2006-2010 is 5 years, then
%each file has 10 years of 12 months, 
mos = 5*12 + 10*12*(length(tstart) - 1);
mod_time = nan(mos,1);
tbnds = nan(2,mos);
tbtm = nan(360,180,mos);

mstart = [1 (1+5*12):120:mos];
mend = 5*12:120:mos;

%%
for t=1:length(tstart)
ncid = netcdf.open([fpath 'gfdl-esm2m_rcp85_to_zb_monthly_',num2str(tstart(t)),'-',...
    num2str(tend(t)),'.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%NaNs on land cells
to(to>=1e19)=nan;
tb = to;

%concatenate
mod_time(mstart(t):mend(t)) = time;
% tbnds(:,mstart(t):mend(t)) = time_bnds;
tbtm(:,:,mstart(t):mend(t)) = tb;

clear to time time_bnds 

end

save([fpath 'gfdl-esm2m_rcp85_tbtm_monthly_',num2str(tstart(1)),'-',...
    num2str(tend(end)),'.mat'],'tbtm','mod_time');
save([cpath 'gfdl-esm2m_rcp85_tbtm_monthly_',num2str(tstart(1)),'-',...
    num2str(tend(end)),'.mat'],'tbtm','mod_time');
