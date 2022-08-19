% Resave daily interp as non-struct

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/Phase3/QuarterDeg/';

%%
yrs = 1961:2010;
for y = 1%:nyrs
    YR = yrs(y);
    
    clear ESM
    load([fpath 'Data_gfdl_mom6_cobalt2_ctrlclim_15arcmin_daily_',num2str(YR),'.mat'],'ESM');
    
    Tp = ESM.Tp;
    Tb = ESM.Tb;
    Zm = ESM.Zm;
    det = ESM.det;
    
    % save
    save([fpath 'gfdl_mom6_cobalt2_ctrlclim_15arcmin_daily_',num2str(YR),'.mat'],...
        'Tp','Tb','Zm','det','-v7.3');
    
end

%%
for y = 1%:nyrs
    YR = yrs(y);
    
    clear ESM
    load([fpath 'Data_gfdl_mom6_cobalt2_obsclim_15arcmin_daily_',num2str(YR),'.mat'],'ESM');
    
    Tp = ESM.Tp;
    Tb = ESM.Tb;
    Zm = ESM.Zm;
    det = ESM.det;
    
    % save
    save([fpath 'gfdl_mom6_cobalt2_obsclim_15arcmin_daily_',num2str(YR),'.mat'],...
        'Tp','Tb','Zm','det','-v7.3');
    
end

%%
tic
ESM=matfile(['/Volumes/MIP/Fish-MIP/Phase3/QuarterDeg/',...
            'gfdl_mom6_cobalt2_obsclim_15arcmin_daily_',ti,'.mat']);
ENVR = get_ESM_empHPl(ESM,GRD,param,DY);
toc

%%
tic
load(['/Volumes/MIP/Fish-MIP/Phase3/QuarterDeg/',...
            'Data_gfdl_mom6_cobalt2_obsclim_15arcmin_daily_',ti,'.mat'],'ESM');
toc
