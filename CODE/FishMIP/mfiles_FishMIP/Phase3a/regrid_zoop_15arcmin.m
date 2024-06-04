% regrid 1/4 files form Xiao to regular 1/4 grid
% to match up with ISIMIP vars for tp, tb, det, etc. 
% Xiao: 1440x1080, ISIMIP: 1440x720

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
cpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';


load([fpath 'gfdl-mom6-cobalt2_15arcmin_mdz_lgz_100m_global_monthly_1961_2010.mat'])
load([fpath 'gfdl-mom6_cobalt2_15arcmin_HPloss_mdz_lgz_100m_month_1961_2010.mat'])