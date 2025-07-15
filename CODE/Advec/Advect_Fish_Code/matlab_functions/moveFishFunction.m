function [ Flux ] = moveFishFunction( Predator, Flux, Prey, apparent_current_full, idx, dt, dx_m, dy_m )
%% Calculates the number of fish based on current and food at cell idx.
% Basically, extracts a single location and calculates the movement
%
% INPUTS:
%   Predator              : [m x n] array of fish biomass  [g/m^3]
%   Flux                  : [m x n] array of flux of biomass due to movement
%   Prey                  : [m x n] array of available food [g/m^3]
%   apparent_current_full : [m x n x 4] apparent current (swim speed taken
%                            into account ) [m/s]
%   idx                   : index of the grid cell [m_i, n_i]
%   dt                    : timestep [s]
%   dx_m, dy_m            : grid spacing in x and y directions [m]
%
% OUTPUTS:
%   Flux                  : updated flux of fish array for grid cell idx.
%
% AUTHOR: [JARED BRZENSKI
% DATE  : 30-06-2025
% -------------------------------------------------------------------------
%% The core of the move fish code.
i = idx(1);
j = idx(2);

[food_dir] = MoreFoodFourWay( Prey, idx );

apparent_current_single = squeeze(apparent_current_full(i,j,:));

Flux = semiLagrangianFish( Predator, Flux, [i,j], food_dir, ...
                         apparent_current_single, ...
                         dt, dx_m, dy_m);
end