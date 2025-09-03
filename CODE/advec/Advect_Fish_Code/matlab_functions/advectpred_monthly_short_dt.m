% Input data and params needed in advection-diffusion scheme

clear 
close all

% Add matlab functions to path
addpath('matlab_functions');

% Velocities
vpath = '../data/';
load([vpath 'feb152013_run25_ocean.198801-200712_uh200_vh200.mat'],'u200','v200', 'geolat_t', 'geolon_t');

% Grid
load([vpath 'Data_hindcast_grid_cp2D.mat'])

% Neighbors
load([ 'all_neighbors_2_360x200.mat'])

% Alt_Velocities
load('Vel100_esm2m_core_daily_1988.mat');

lon = geolon_t;
lat = geolat_t;

v = VideoWriter('fish_anim_swim_w_negs_1980_1.avi');
v.FrameRate = 10;  % Frames per second
open(v);
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
bio = ones(ni,nj);
%bio = 11*rand(ni,nj) - 1;
%Global
%bio = 10*ones(ni,nj);   %Global
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
%prey = zeros(ni,nj);
%prey = 10*rand(ni,nj);
prey = 0 + ( 10 - 0 ) * cosd(abs(lat) / 2);

%prey = 100*ones(ni,nj);   %Global
%prey(220:240,:) = 10.0; prey(121:141,195:200) = 10.0; %Atl-Arctic
%prey(:,84:109) = 1.0e1;     %seed equator
%prey(220:240,:) = 1.0e1;    %seed Atl
%prey(59:79,:) = 1.0e1;      %seed Pac
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
tstep = 24 * 60 * 60 / 4; %time step in seconds

% Files to save
cname='Atl_even_dt1d_velMO2_b100';
biov = zeros(NX,DAYS*YEARS);
preyv = prey(ID);

fish_speed = 1.0; %(m/s)

alt_current = zeros(ni,nj,2);

%% call advec-diff
M=0;
n=0;
fig = figure('Units', 'pixels', 'Position', [100 100 560 420]);  % match required size

for Y=1:YEARS
    for mo = 1:length(Mos)
        M = M+1;
        alt_current_u = ESM.U(:,M);
        alt_current_v = ESM.V(:,M);
        temp = zeros(ni,nj);
        temp(ID) = alt_current_u;
        alt_current(:,:,1) = temp;
        temp = 0.*temp;
        temp(ID) = alt_current_v;
        alt_current(:,:,2) = temp;
        %current(:,:,1) = u200(:,:,M); 
        %current(:,:,2) = v200(:,:,M);
        for DAY = 1:Mos(mo) 
            step_bio = bio(ID);
            n=n+1;
            [num2str(mo) ',' num2str(DAY)];
            for tt= 1:4
                bio = AdvectPredator( bio,...
                                  prey, ...
                                  alt_current, ...
                                  tstep, ...
                                  GRD.dxtn, ...
                                  GRD.dyte, ...
                                  neighborhood, ...
                                  fish_speed, ...
                                  GRD.mask, ...
                                  GRD.area, ...
                                  nj, ...% m from above
                                  ni);   % n from above
            end

            biov(:,n) = bio(ID); 
            STEP_sum = sum(sum( step_bio .* GRD.area(ID) ) );
            FINAL_sum = sum(sum( bio(ID) .* GRD.area(ID) ) ); 
            %fprintf('day: %2d, mnth: %2d, step_diff: %5.6e, overall_diff: %5.6e \n',DAY, mo, OG_sum-STEP_sum, OG_sum-FINAL_sum);
            %fprintf('230,98 = %4.5e\n', bio(232,95));
            %fprintf('Max velocity: %4.5e\n', max(max(max(alt_current))));
                %% Plot Fish Concentration
            %set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
            ax1 = gca;
            img1 = pcolor(lon, lat, bio);
            img1.FaceColor='flat';
            img1.EdgeColor='none';
            cmap1 = cmocean('dens');
            colormap(ax1,cmap1); cb1 = colorbar;
            ylabel(cb1, '$gww/m^3$','FontSize',14,'interpreter','latex');
            set(cb1, 'TickLabelInterpreter','latex');
            axis tight
            pbaspect([2 1 1])
            xlabel(ax1,'x [km]', 'interpreter','latex');
            ylabel(ax1,'y [km]', 'interpreter','latex');
            title(['Fish Con on Day: ', num2str(n)]);
            clim([0 10]);
            drawnow
            %frame = getframe(fig);
            %writeVideo(v, frame);
        end
    end
end
close(v);
%% Save
spath = 'data/';
save([spath 'AdvectPred2_01_speed' cname '.mat'],'biov','preyv','GRD');

