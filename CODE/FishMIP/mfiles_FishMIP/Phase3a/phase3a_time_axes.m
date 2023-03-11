% Create "days since" start date vectors for Phase 3a sims

clear 
close all

load('phase3a_time_axis.mat')

%% create datenum to extend backwards 10 yrs
dnum_ymd = datenum(Y,M,D);

test=1831:1840;
Yspin=repelem(test,12);
test3=1:12;
nyr = length(Yspin)/12;
Mspin=repmat(test3,1,nyr);
Dspin = ones(1,length(Yspin));

dnum_spin = datenum(Yspin',Mspin',Dspin');

%% Historic 1961-2010
% time is supposed to be days since Jan 1, 1901
% Cami gave us this, just need to take subset
hid = find(Y>=1961);
hist_time = dsince(hid);
hist_units = 'days since Jan-01-1901';

%% Spinup 1831-1840
% just want last 10 yrs
% time can be days since Jan 1, 1831
% use 1st 10 yrs since 1901
sid = find(Y>=1901);
spin_time = dsince(sid(1:length(Yspin)));
spin_units = 'days since Jan-01-1831';

%% Transition (PI) 1841-1960
% time can be days since Jan 1, 1841
% shift the csv that Cami sent
pid = find(Y==1841);
dtime = dsince(pid(1));
pi_dsince = dsince + abs(dtime);
pid2 = find(Y>=1841 & Y<1961);
pi_time = pi_dsince(pid2);
pi_units = 'days since Jan-01-1841';

%%
save('FishMIP_phase3a_exper_times.mat','spin_time','pi_time','hist_time',...
    'spin_units','pi_units','hist_units');
