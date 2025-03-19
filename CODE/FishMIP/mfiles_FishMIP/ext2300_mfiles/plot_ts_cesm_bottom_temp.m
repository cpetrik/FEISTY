% SSP 534-over 2101-2299
% plot theta-bot because looks off from tob in other scenarios

clear 
close all

%%
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp534over/';

load([fpath 'cesm2_ssp534-over_temp_btm_monthly_2040_2299.mat']);
load([fpath 'cesm2-waccm_r1i1p1f1_ssp534-over_deptho_60arcmin_global_fx.mat'],'deptho')

temp_btm(temp_btm > 1.0e19) = nan;

ssp534_Tb = temp_btm;
ssp534_dep = deptho;

clear deptho temp_btm

%% 
spath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp585/';

load([spath 'cesm2_ssp585_temp_btm_monthly_2015_2299.mat']);
load([spath 'cesm2-waccm_r1i1p1f1_ssp585_deptho_60arcmin_global_fx.mat'],'deptho')

temp_btm(temp_btm > 1.0e19) = nan;

ssp585_Tb = temp_btm;
ssp585_dep = deptho;

clear deptho temp_btm

%%
mos = length(time);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

yrs = floor(yr(1)):floor(yr(end));
Tdays=1:365;

%%
