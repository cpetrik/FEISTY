% Input data and params needed in advection-diffusion scheme

clear
close all

% Add matlab functions to path
addpath('matlab_functions');

% Velocities
%vpath = 'data/';
vpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CORE-forced/';
load([vpath,'Vel100_esm2m_core_daily_1988.mat'],'ESM');

% Grid
load([vpath 'Data_hindcast_grid_cp2D.mat'])
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');

% Neighbors
load(['all_neighbors_2_360x200.mat'])

lon = geolon_t;
lat = geolat_t;

%% Initial fish and prey biomass

load('/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/CORE/Last_mo_Spinup1988_no_move_All_fish03_Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100.mat')

load([vpath,'Data_ocean_cobalt_daily_1988.mat'],'COBALT');

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
bio = nan*ones(ni,nj);
bio(ID) = Med_f.bio;
bio = bio .* GRD.mask;

OG_sum = sum(sum( bio(ID) .* GRD.area(ID) ) );

%% define prey
prey = zeros(ni,nj);
prey(ID) = COBALT.Zl(:,1);
prey = prey .* GRD.mask;

%%
figure
pcolor(geolat_t,geolon_t,bio); shading flat

figure
pcolor(geolat_t,geolon_t,prey); shading flat

%% define time
YEARS = 1;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
Mos = repmat(MNTH,1,YEARS);
tstep = 6 * 60 * 60; %time step in seconds

% Files to save
cname='Global_ICs_MF_LZ_dt6h_velDAY_swim01';
biov = zeros(NX,DAYS*YEARS);
preyv = prey(ID);

fish_speed = 0.10; %(m/s)

%% call advec-diff
M=0;
n=0;
for Y=1:YEARS
    for DAY = 1:365
        Utemp = nan*ones(ni,nj);
        Vtemp = nan*ones(ni,nj);
        Utemp(ID) = ESM.U(:,DAY);
        Vtemp(ID) = ESM.V(:,DAY);

        current(:,:,1) = Utemp;
        current(:,:,2) = Vtemp;
        
        %prey = prey.*rand(ni,nj);

        step_bio = bio(ID);
        n=n+1;
        %[num2str(mo) ',' num2str(DAY)];
       % num2str(DAY)
        bio = AdvectPredator( bio,...
            prey, ...
            current, ...
            tstep, ...
            GRD.dxtn, ...
            GRD.dyte, ...
            neighborhood, ...
            fish_speed, ...
            GRD.mask, ...
            GRD.area, ...
            nj, ...
            ni);

        biov(:,n) = bio(ID);
        %STEP_sum = sum(sum( step_bio .* GRD.area(ID) ) );
        %FINAL_sum = sum(sum( bio(ID) .* GRD.area(ID) ) );
        %fprintf('day: %2d, mnth: %2d, step_diff: %5.6e, overall_diff: %5.6e \n',DAY, mo, OG_sum-STEP_sum, OG_sum-FINAL_sum);
        %fprintf('230,98 = %4.5e\n', bio(232,95));

    end
end

%% Save
spath = '/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/CORE/';
save([spath 'AdvectPred_' cname '.mat'],'biov','preyv','GRD');

