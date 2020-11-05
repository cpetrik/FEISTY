% Compare timeseries of GFDL and IPSL temp and zoop
% Just use surface temp and zmeso vint as indicators

clear all
close all

ipath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/';
gpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/';

%%
load('zoop_ts.mat')

%%
figure(1)
plot(ghts,'b'); hold on;
plot(ihts,'r'); hold on;

figure(2)
plot(ghzm,'b'); hold on;
plot(ihzm,'r'); hold on;

%%
mean(ghts)
mean(ihts)
mean(ghzm)
mean(ihzm)

%%
figure(3)
subplot(2,2,1)
histogram(ghts)
subplot(2,2,2)
histogram(ihts)
subplot(2,2,3)
histogram(ghzm)
subplot(2,2,4)
histogram(ihzm)

%%
figure
plot(ipts,'r'); hold on;
