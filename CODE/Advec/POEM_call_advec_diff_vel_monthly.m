% Input data and params needed in advection-diffusion scheme

clear all
close all

% Velocities
% vpath = '/Volumes/GFDL/GCM_DATA/CORE-forced/';
% load([vpath 'feb152013_run25_ocean.198801-200712_uh200_vh200.mat'],'u200','v200');
vpath = '/Volumes/GFDL/GCM_DATA/ESM2M_hist/';
load([vpath 'ocean.199601-200012_uh200_vh200.mat'],'u200','v200');

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

%% define a patch to advect
bio = zeros(ni,nj);
%Global
bio = 100*ones(ni,nj);   %Global
%bio(220:240,:) = 10.0; bio(121:141,195:200) = 10.0; %Atl-Arctic
%bio(:,84:109) = 1.0e1;     %seed equator
%bio(220:240,:) = 1.0e1;    %seed Atl
%bio(59:79,:) = 1.0e1;      %seed Pac
%bio(5:25,:) = 1.0e1;       %seed Indian W
%bio(340:360,:) = 1.0e1;    %seed Indian E
%bio(:,181:200) = 1.0e1;    %seed Arctic
%bio(:,12:32) = 1.0e1;      %seed Antarctic

bio = bio .* GRD.mask;

% define diffusivity
K = 600.0;

% define time
YEARS = 1;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
Mos = repmat(MNTH,1,YEARS);
tstep = 1; %time step in hours

% Files to save
cname='Global_even_dt1hr_esm2m2000_velMO_b100_area';
biov = zeros(NX,DAYS*YEARS);

%% call advec-diff
%M=0;
%Test 2000 = last year
M=48;
n=0;
for Y=1:YEARS
    for mo = 1:length(Mos)
        M = M+1;
        U = u200(:,:,M); 
        V = v200(:,:,M);
        for DAY = 1:Mos(mo)
            n=n+1;
            [num2str(mo) ',' num2str(DAY)]
            bio = sub_advec_diff_vel(GRD,bio,K,U,V,ni,nj,tstep);
            biov(:,n) = bio(ID);
        end
    end
end

% Save
csvwrite(['/Volumes/GFDL/CSV/advect_tests/Matlab_adv_diff_' cname '.csv'],biov);

