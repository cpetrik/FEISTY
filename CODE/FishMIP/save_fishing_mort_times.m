% Divide fishing mortality into time periods

clear all 
close all

%% 1961-2010
spath = '/Volumes/MIP/Fish-MIP/Phase3/fishing/grid_mortality_guilds/';
cpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/grid_mortality_guilds/';
fpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/fishing_ms_ideas/fishing_effort_impl/grid_mortality_guilds/';
load([fpath 'grid_mortality_all.mat'])
save([spath 'grid_mortality_all.mat'])
save([cpath 'grid_mortality_all.mat'])

%% two time periods
%transition = 1841-1960
%historic = 1961-2010

yt = (yrD<1961);
yh = (yrD>=1961);

fmortDt = fmortD(:,yt);
fmortDh = fmortD(:,yh);
fmortFt = fmortF(:,yt);
fmortFh = fmortF(:,yh);
fmortPt = fmortP(:,yt);
fmortPh = fmortP(:,yh);

yrt = yrD(yt);
yrh = yrD(yh);

%% save indep
fmortD = fmortDt;
fmortF = fmortFt;
fmortP = fmortPt;
yr = yrt;

save([fpath 'grid_mortality_all_1841_1960.mat'],'fmortF','fmortP','fmortD','yr',...
    'LatF','LonF','LatP','LonP','LatD','LonD')
save([spath 'grid_mortality_all_1841_1960.mat'],'fmortF','fmortP','fmortD','yr',...
    'LatF','LonF','LatP','LonP','LatD','LonD')
% save([cpath 'grid_mortality_all_1841_1960.mat'],'fmortF','fmortP','fmortD','yr',...
%     'LatF','LonF','LatP','LonP','LatD','LonD')

clear fmortF fmortD fmortP yr

%%
fmortD = fmortDh;
fmortF = fmortFh;
fmortP = fmortPh;
yr = yrh;

save([fpath 'grid_mortality_all_1961_2010.mat'],'fmortF','fmortP','fmortD','yr',...
    'LatF','LonF','LatP','LonP','LatD','LonD')
save([spath 'grid_mortality_all_1961_2010.mat'],'fmortF','fmortP','fmortD','yr',...
    'LatF','LonF','LatP','LonP','LatD','LonD')
% save([cpath 'grid_mortality_all_1961_2010.mat'],'fmortF','fmortP','fmortD','yr',...
%     'LatF','LonF','LatP','LonP','LatD','LonD')



