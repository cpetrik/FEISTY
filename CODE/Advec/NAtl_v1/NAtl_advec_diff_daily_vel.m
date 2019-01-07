% Advection and diffusion using upwind scheme adapted from MOM code
% Specifically in the N Atl region of N Henschke jellyfish model
% This loads and defines all the constants and parameters you need for
% advection and diffusion
% Then it calls the advection-diffusion or pure advection function
% You can put this in the main body of your code
% Colleen M. Petrik 2/15/17

clear all
close all

% Grid
load('Data_hindcast_grid2D.mat')

%% Natasha's region
iids = [206:295];
jids = [150:177];
% Add one extra cell on each side to make i+1, i-1, j+1, j-1 calcs easy
ni = length(iids)+2;
nj = length(jids)+2;
gis = iids(1)-1;
gjs = jids(1)-1;
gie = iids(end)+1;
gje = jids(end)+1;

%% Define a patch to advect (this would be your biomasses, etc.)
bio = zeros(ni,nj);
bio(2:ni-1,2:nj-1) = 100;   %Put your [90x28] matrix here

bio = bio .* GRD.mask(gis:gie,gjs:gje);

% define diffusivity
K = 600.0;

% define time
YEARS = 1;
DAYS = 365;
tstep = 1; %time step in hours

% Files to save
NX = length(iids)*length(jids);
biov = zeros(NX,DAYS*YEARS);
cname='NAtl_dt1hr_vel50_daily_b100';

%% Do advec-diff
% This can be moved into your daily time step loop

n=0;
for Y=1:YEARS
    yr = num2str(Y+1948-1);
    % Velocities
    load(['Vel50_feb152013_run25_ocean_' yr '.mat'],'u','v');
    for DAY = 1:DAYS
        DAY
        n=n+1;
        U = u(:,:,DAY); 
        V = v(:,:,DAY);
        bio = sub_NAtl_advec_diff_vel(GRD,bio,K,U,V,iids,jids,tstep);
        bio2 = bio(2:ni-1,2:nj-1);
        biov(:,n) = reshape(bio2,length(iids)*length(jids),1);
    end
end

% Save biomass each day (for testing)
csvwrite(['Adv_diff_' cname '.csv'],biov);

%% Do advec only
% % This can be moved into your daily time step loop
% 
% n=0;
% for Y=1:YEARS
%     yr = num2str(Y+1948-1);
%     % Velocities
%     load(['Vel50_feb152013_run25_ocean_' yr '.mat'],'u','v');
%     for DAY = 1:DAYS
%         DAY
%         n=n+1;
%         U = u(:,:,DAY); 
%         V = v(:,:,DAY);
%         bio = sub_NAtl_advec_vel(GRD,bio,U,V,iids,jids,tstep);
%         bio2 = bio(2:ni-1,2:nj-1);
%         biov(:,n) = reshape(bio2,length(iids)*length(jids),1);
%     end
% end
% 
% % Save
% csvwrite(['Adv_' cname '.csv'],biov);

