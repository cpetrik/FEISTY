clear all

%%
load(['/Volumes/MIP/Fish-MIP/Phase3/OneDeg/',...
    'Data_gfdl_mom6_cobalt2_ctrlclim_onedeg_daily_1961.mat']);
Tp=ESM.Tp;
Tb=ESM.Tb;
det=ESM.det;
Zm=ESM.Zm;

save(['/Volumes/MIP/Fish-MIP/Phase3/OneDeg/',...
    'test_Data_gfdl_mom6_cobalt2_ctrlclim_onedeg_daily_1961.mat'],...
    'Tp','Tb','Zm','det','-v7.3');

%%
clear all
egESM=matfile(['/Volumes/MIP/Fish-MIP/Phase3/OneDeg/',...
    'test_Data_gfdl_mom6_cobalt2_ctrlclim_onedeg_daily_1961.mat']);

%%
DY=1;
ENVR.Tp(:,1)  = egESM.Tp(:,DY);
ENVR.Tb(:,1)  = egESM.Tb(:,DY);
ENVR.Zm(:,1)  = egESM.Zm(:,DY);
ENVR.det(:,1) = egESM.det(:,DY);
% HP loss is an empirical fitted fn of biomass and temp
ENVR.dZm(:,1) = 10 .^ (-2.925 + 1.964.*log10(ENVR.Zm(:,1)+eps) + 1.958e-2.*ENVR.Tp(:,1));

%%
%param.ID = length()
ENVR.fZl(:,1) = zeros(param.NX,1);
ENVR.fB(:,1)  = zeros(param.NX,1);




