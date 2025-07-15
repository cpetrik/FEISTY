%% Simple test case for AdvectPredator
%  Fish follow food
%
% Test Case 1: Linear advection for 3 days.
% No consumption, no speed adjustment, just
% if Conc_n > Conc_m, swim towards Cell_n at max_speed
%
% A lot of this is to set up a domain, and dt, fish pop, etc.
% The actual function call starts at line 101.
%
% Functions NOT needed by the AdvectPredator.m:
%     makeFood.m
%     zeroCornerMatrix.m
% 
% BY: Jared Brzenski
%
%
%%
clc 
close all
clearvars

% Add matlab_functions to path
addpath(genpath('matlab_functions'));

%% Time / draw
days = 3;             % sim time
num_hours = 2.25;     % dt
iters = ceil(days * 24 / num_hours);             

%% System Parameters
% Fish Swimming Speed
% Here BATS Tuna Speed 86-260km/day ~ 1m/s to 3m/s
% Biomass is approximately 1gww/m^3 ~ 
% Tuna consume ~5% of their mass per day
fish_speed = 1.0;       % speed in meters/second
%
fish_init = 1.0;        % fish biomass in g / m^3
food_init = 1000;       % food biomass in g / m^3
%
% Currents
u_max = 0.;             % current speed in meters/second
v_max = 0.0;
% Time Discretization
dt   = num_hours/24;    % dt in days
dt_s = dt*24*60*60;     % dt in seconds
day  = 24*60*60;        % seconds in a day ( for dimensionalizing )
%
epsilon = 10e-15;      % tiny value, 
%
%% Grid setup
% Number of CELLS in each direction
m = 100;
n = 50;
% Domain's start and end locations ( distances )
a = 0;
b = 1000;
c = -250.0;
d = 250.0;
% Spatial step sizes
dx = (b-a)/(m-1);         % kilometers
dy = (d-c)/(n-1);         % kilometers
dx_m = dx * 1000.0;   % dx in meters
dy_m = dy * 1000.0;   % dy in meters
% Grid
[X, Y] = meshgrid( a : dx : b, c : dy : d );
%
%% Initial Density Generation for Fish
Fish = fish_init * ones(n, m);
Fish = zeroCornerMatrix(Fish);
Fish_init = Fish;
%
x0_fish = 300; % km
Fish = makeFood('humpfish', Fish_init, X, Y, x0_fish, 0, 1, 0);

%% Initial Food Generation
Food = food_init * ones(n, m);
Food = makeFood('gradient', Food, X, Y, 0, 0, 1, 0);
Food = fliplr(Food); % Flip gadient left to right
Food = Food .* food_init;

%% Island Masking
mask = ones(size(Food));
Fish = mask .* Fish;
Food = mask .* Food; 

%% Initial Current Field
current = zeros(n,m,2);

%% Convert to day from seconds
% Convert to m per day
current = current .* day;
% Fish speed to m per day
fish_speed = fish_speed .* day;

%% Mask Values
current(:,:,1) = squeeze(current(:,:,1)) .* mask;
current(:,:,2) = squeeze(current(:,:,2)) .* mask;
apparent_current = current;

Food = Food .* mask;
Fish = Fish .* mask;

%% Iteration Starts
for it = 1 : iters
    [ Fish ] = AdvectPredator( Fish,...
                              Food, ...
                              current, ...
                              dt, ...
                              dx_m, ...
                              dy_m, ...
                              fish_speed, ...
                              mask, ...
                              m, ...
                              n);

    fprintf('Iteration: %5d/%5d\n', it, iters);
    
    %% Plot Fish Concentration
    figure(23);
    ax1 = gca;
    img1 = pcolor(X, Y, Fish);
    img1.FaceColor='interp';
    img1.EdgeColor='none';
    cmap1 = cmocean('dens');
    colormap(ax1,cmap1); cb1 = colorbar;
    ylabel(cb1, '$gww/m^3$','FontSize',14,'interpreter','latex');
    set(cb1, 'TickLabelInterpreter','latex');
    axis tight
    pbaspect([2 1 1])
    xlabel(ax1,'x [km]', 'interpreter','latex');
    ylabel(ax1,'y [km]', 'interpreter','latex');
    title(ax1,'Fish Concentration in $gww/m^3$', 'interpreter','latex');
    %clim([0. 0.1]);
    caxis([0 0.1]);
    %
end

