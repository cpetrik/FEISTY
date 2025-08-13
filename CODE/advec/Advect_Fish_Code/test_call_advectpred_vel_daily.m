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
bio = 10*ones(ni,nj);   %Global
%bio(220:240,:) = 10.0; bio(121:141,195:200) = 10.0; %Atl-Arctic
%bio(:,84:109) = 1.0e1;     %seed equator
%bio(220:240,:) = 1.0e1;    %seed Atl
%bio(59:79,:) = 1.0e1;      %seed Pac
%bio(5:25,:) = 1.0e1;       %seed Indian W
%bio(340:360,:) = 1.0e1;    %seed Indian E
%bio(:,181:200) = 1.0e1;    %seed Arctic
%bio(:,12:32) = 1.0e1;      %seed Antarctic

bio = bio .* GRD.mask;

OG_sum = sum(sum( bio(ID) .* GRD.area(ID) ) );

%% define prey
prey = zeros(ni,nj);
%prey = 100*ones(ni,nj);   %Global
%prey(220:240,:) = 10.0; prey(121:141,195:200) = 10.0; %Atl-Arctic
%prey(:,84:109) = 1.0e1;     %seed equator
prey(220:240,:) = 1.0e1;    %seed Atl
prey(59:79,:) = 1.0e1;      %seed Pac
%prey(5:25,:) = 1.0e1;       %seed Indian W
%prey(340:360,:) = 1.0e1;    %seed Indian E
%prey(:,181:200) = 1.0e1;    %seed Arctic
%prey(:,12:32) = 1.0e1;      %seed Antarctic

prey = prey .* GRD.mask;

%% define time
YEARS = 1;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
Mos = repmat(MNTH,1,YEARS);
tstep = 24 * 60 * 60; %time step in seconds

% Files to save
cname='Global_evenfish_evenprey_dt1d_velDAY_b100_swim01';
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

