% Plot BOATS & FEISTY against obs catches
% In US LMEs for MAPP proposal

clear
close all

%% Fishing data
%fpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

load([fpath 'FishMIP_Phase3a_LME_Catch_1948-2015_interann_var.mat'],...
    'lme_a_mean','lme_a_ts','units');

%% BOATS output
bpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/BOATS/';

load([bpath 'boats_gfdl-mom6-cobalt2_obsclim_histsoc_tc_LME_annual_1961_2010.mat'],...
    'lme_mcatch');

blme_mcatch = lme_mcatch * 1e-6; % g -> MT
clear lme_mcatch

%% FEISTY output
ppath='/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100/QuarterDeg/';

load([ppath 'feisty_gfdl-mom6-cobalt2_obsclim_histsoc_tc_LME_annual_1961_2010.mat'],...
    'lme_mcatch');

plme_mcatch = lme_mcatch * 1e-6; % g -> MT
clear lme_mcatch

%% Align times
%Catch 1948-2015
%MEMs 1961-2010

cyr = 1948:2015;
myr = 1961:2010;

[yr,cid] = intersect(cyr,myr);

lme_a_ts = lme_a_ts(:,cid);

%% figs
pp = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/FishMIP/Phase3a/';
if (~isfolder(pp))
    mkdir(pp)
end

%% NPac and NEAtl LMEs
% lid = [54,1,2,65,3,6,7]; 
% lname = {'CHK','EBS','GAK','AI','CCE','SE','NE'};
% CHK in arctic grid
lid = [1,2,65,3,6,7]; 
lname = {'EBS','GAK','AI','CCE','SE','NE'};

f1 = figure('Units','inches','Position',[3 3 9 6]);
for i=1:length(lid)
    subplot(3,2,i)
    plot(yr,1e-6*lme_a_ts(lid(i),:),'k'); hold on;
    plot(yr,1e-6*blme_mcatch(lid(i),:),'b'); hold on;
    plot(yr,1e-6*plme_mcatch(lid(i),:),'r'); hold on;
    ylabel('mil MT')
    if (i>3)
        xlabel('Year')
    end
    title(lname{i})
    if (i==6)
        legend('obs','BOATS','FEISTY')
    end
end
print('-dpng',[pp 'ts_USLMEs_tc_boats_feisty_obs.png'])

%% Plot as anomalies from mean
aclme_mcatch = lme_a_ts - mean(lme_a_ts,2);
ablme_mcatch = blme_mcatch - mean(blme_mcatch,2);
aplme_mcatch = plme_mcatch - mean(plme_mcatch,2);

f2 = figure('Units','inches','Position',[3 3 9 6]);
for i=1:length(lid)
    subplot(3,2,i)
    plot(yr,1e-6*aclme_mcatch(lid(i),:),'k'); hold on;
    plot(yr,1e-6*ablme_mcatch(lid(i),:),'b'); hold on;
    plot(yr,1e-6*aplme_mcatch(lid(i),:),'r'); hold on;
    title(lname{i})
    ylabel('mil MT')
    if (i>3)
        xlabel('Year')
    end
    if (i==4)
        legend('obs','BOATS','FEISTY')
        legend('location','northwest')
    end
end
print('-dpng',[pp 'ts_anom_USLMEs_tc_boats_feisty_obs.png'])
