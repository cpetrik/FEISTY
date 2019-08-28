% Read Fish-MIP CESM Historic netcdfs
% Integrate over top 100 m for small zooplankton

clear all
close all

fpath='/Volumes/FEISTY/Fish-MIP/GFDL/Hist/';
cpath='/Volumes/GFDL/Fish-MIP/GFDL/Hist/';

tstart = 1861:10:2001;
tend = 1870:10:2005;
tend = [tend 2005];

%each file has 10 years, except last has 5
mos = 10*(length(tstart)-1) + 5;
mod_time = nan(mos,1);
tbnds = nan(2,mos);
nsmz_100 = nan(360,180,mos);

mstart = 1:10:mos;
mend = 10:10:mos;
mend = [mend mos];

%%
for t=1:length(tstart)
ncid = netcdf.open([fpath 'gfdl-esm2m_hist_szoo_zall_annual_',num2str(tstart(t)),'-',...
    num2str(tend(t)),'.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%NaNs on land cells
szoo(szoo>=9e+36)=nan;
%top 100 m
sz100 = szoo(:,:,1:10,:);
%integrated over top 100m
isz100 = squeeze(nansum(sz100,3));

%concatenate
mod_time(mstart(t):mend(t)) = time;
tbnds(:,mstart(t):mend(t)) = time_bnds;
nsmz_100(:,:,mstart(t):mend(t)) = isz100;

clear szoo time time_bnds 

end

save([fpath 'gfdl-esm2m_hist_szoo_100_monthly_',num2str(tstart(1)),'-',...
    num2str(tend(end)),'.mat'],'nsmz_100','mod_time','tbnds');
save([cpath 'gfdl-esm2m_hist_szoo_100_monthly_',num2str(tstart(1)),'-',...
    num2str(tend(end)),'.mat'],'nsmz_100','mod_time','tbnds');


