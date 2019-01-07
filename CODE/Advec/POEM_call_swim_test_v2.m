% Input data and params needed in advection-diffusion scheme

clear all
close all

tic

% Transports path
% vpath = '/Volumes/GFDL/GCM_DATA/CORE-forced/';
vpath = '/Volumes/GFDL/POEM_JLD/esm2m_hist/';

% Grid
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Data_hindcast_grid_cp2D.mat')

%% number of water cells
ID = find(GRD.mask==1);
NX = length(ID);

% grid size
[ni,nj] = size(GRD.mask);
isd = 1;
jsd = 1;
ied = ni;
jed = nj;

% define a patch to advect
bio = zeros(ni,nj);
%Global
bio = 100*ones(ni,nj);   %Global
%bio(220:240,:) = 100.0; bio(121:141,195:200) = 100.0; %Atl-Arctic
%bio(:,84:109) = 1.0e2;     %seed equator
%bio(220:240,:) = 1.0e2;    %seed Atl
%bio(59:79,:) = 1.0e2;      %seed Pac
%bio(5:25,:) = 1.0e2;       %seed Indian W
%bio(340:360,:) = 1.0e2;    %seed Indian E
%bio(:,181:200) = 1.0e2;    %seed Arctic
%bio(:,12:32) = 1.0e2;      %seed Antarctic
%bio(206:295,150:177) = 1.0e2; %w/i Natasha Atl

bio = bio .* GRD.mask;

% define time
YEARS = 1;
DAYS = 365;
tstep = 1; %time step in hours

%% Swimming behavior

% load('/Volumes/GFDL/GCM_DATA/CORE-forced/ocean_cobalt.feb15_run25.1988-2007_phyt_zoop.mat',...
%     'nmdz_avg200_88','nlgz_avg200_88')
% load('/Volumes/GFDL/GCM_DATA/CORE-forced/ocean.186101-200512.temp_100_avg.mat',...
%     'TEMP_100')
% zoop = nmdz_avg200_88(:,:,end-11:end) + nlgz_avg200_88(:,:,end-11:end);
% T100 = TEMP_100(:,:,end-11:end) - 273.15;
% mid = find(T100 == min(T100(:)));
% T100(mid) = NaN;
% zid = find(zoop == min(zoop(:)));
% zoop(zid) = NaN;
% jdmo = 15:30:365;
% zoopi = zeros(ni,nj,365);
% ti = zoopi;
% for i=1:ni
%     for j=1:nj
%         zoopi(i,j,:) = interp1(jdmo,squeeze(zoop(i,j,:)),1:365,'linear','extrap');
%         ti(i,j,:) = interp1(jdmo,squeeze(T100(i,j,:)),1:365,'linear','extrap');
%     end
% end
% clear TEMP_100 nmdz_avg200_88 nlgz_avg200_88
load('/Volumes/GFDL/GCM_DATA/CORE-forced/swim_test_zoop_temp.mat');
zoopi(isnan(zoopi)) = 0.0;
ti(isnan(ti)) = 0.0;

% swimming speed 
L_m = 200.0; % medium
L_l = 1.0e3; % large
w = exp(0.063*(ti-15.0)) * 0.5*L_m*1e-3;  %medium
%w = exp(0.063*(ti-15.0)) * 0.5*L_l*1e-3;  %large

% define value to maximize
%nu = GRD.ht;           %go to deep
%nu = -1.0 * GRD.ht;    %go to shallow
%nu = zoopi(:,:,75);    %go to food
nu = ti(:,:,75);       %go to high temp

% files to save
cname='Global_even_dt1hr_esm2m2000_vel_b100_area_temp_const';
biov = zeros(NX,DAYS*YEARS);

%% do swim + advec-diff
% % define diffusivity
% K = 600.0;
% n=0;
% for Y=1:YEARS
%     yr = num2str(Y+1988-1);
%     % Velocities
%     %load([vpath 'Vel200_ESM2Mhist_' num2str(yr) '.mat'],'u','v');
%     load([vpath 'Vel200_ESM2Mhist_2000.mat'],'u','v');
%     for DAY = 1:DAYS
%         DAY
%         n=n+1;
%         U = u(:,:,DAY); 
%         V = v(:,:,DAY);
%         Q = w(:,:,DAY) .* GRD.mask;
%         %nu = zoopi(:,:,DAY);     %go to food
%         %nu = ti(:,:,DAY);       %go to high temp
%         bio = sub_swimgrad_advec_diff_vel(GRD,bio,K,U,V,nu,Q,ni,nj,tstep);
%         biov(:,n) = bio(ID);
%     end
% end
% 
% % Save
% csvwrite(['/Volumes/GFDL/CSV/advect_tests/Matlab_swimgrad_adv_diff_' cname '.csv'],biov);

%% do swim + advec only
% % define diffusivity
% K = 0.0;
% n=0;
% for Y=1:YEARS
%     yr = num2str(Y+1988-1);
%     % Velocities
%     %load([vpath 'Vel200_ESM2Mhist_' num2str(yr) '.mat'],'u','v');
%     load([vpath 'Vel200_ESM2Mhist_2000.mat'],'u','v');
%     for DAY = 1:DAYS
%         DAY
%         n=n+1;
%         U = u(:,:,DAY); 
%         V = v(:,:,DAY);
%         Q = w(:,:,DAY) .* GRD.mask;
%         %nu = zoopi(:,:,DAY);     %go to food
%         %nu = ti(:,:,DAY);       %go to high temp
%         bio = sub_swimgrad_advec_diff_vel(GRD,bio,K,U,V,nu,Q,ni,nj,tstep);
%         biov(:,n) = bio(ID);
%     end
% end
% 
% % Save
% csvwrite(['/Volumes/GFDL/CSV/advect_tests/Matlab_swimgrad_adv_' cname '.csv'],biov);

%% do swim only
%define diffusivity
K = 0.0;
U = zeros(ni,nj);
V = zeros(ni,nj);
n=0;
for Y=1:YEARS
    yr = num2str(Y+1988-1);
    for DAY = 1:DAYS
        DAY
        n=n+1;
        Q = w(:,:,DAY) .* GRD.mask;
        %nu = zoopi(:,:,DAY);     %go to food
        %nu = ti(:,:,DAY);       %go to high temp
        bio = sub_swimgrad_advec_diff_vel(GRD,bio,K,U,V,nu,Q,ni,nj,tstep);
        biov(:,n) = bio(ID);
    end
end

% Save
csvwrite(['/Volumes/GFDL/CSV/advect_tests/Matlab_swimgrad_' cname '.csv'],biov);



toc


