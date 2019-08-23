% Read Fish-MIP PreIndustrial netcdfs
% Concatenate column-integrated zooplankton

clear all
close all

fpath='/Volumes/FEISTY/Fish-MIP/GFDL/PreIndust/';

tstart = 1861:10:2091;
tend = 1870:10:2100;

%each file has 10 years 
mos = 10*(length(tstart));
mod_time = nan(mos,1);
tbnds = nan(2,mos);
nz_int = nan(360,180,mos);

mstart = 1:10:mos;
mend = 10:10:mos;

%%
for t=1:length(tstart)
ncid = netcdf.open([fpath 'gfdl-esm2m_pi_zoo_zint_annual_',num2str(tstart(t)),'-',...
    num2str(tend(t)),'.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%NaNs on land cells
zoo(zoo>=1e+19)=nan;
zint = zoo;

% save([fpath 'gfdl-esm2m_pi_zoo_zint_annual_',num2str(tstart(t)),'-',...
%     num2str(tend(t)),'.mat'],'zint','time','time_bnds');

%concatenate
mod_time(mstart(t):mend(t)) = time;
tbnds(:,mstart(t):mend(t)) = time_bnds;
nz_int(:,:,mstart(t):mend(t)) = zint;

clear zoo zint time time_bnds 

end

save([fpath 'gfdl-esm2m_pi_zoo_zint_annual_',num2str(tstart(1)),'-',...
    num2str(tend(end)),'.mat'],'nz_int','mod_time','tbnds');


