function[ Predator ] = AdvectPredator( Predator,...
                                        Prey, ...
                                        current, ...
                                        dt, ...
                                        dx_m, ...
                                        dy_m, ...
                                        fish_speed, ...
                                        mask, ...
                                        m, ...
                                        n)
%% AdvectPredator.m
% -------------------------------------------------------------------------
% Simulates fish advection and feeding in a 2D field using a 
% semi-Lagrangian approach.
%
% This function updates the fish biomass field (`Fish`) based on currents,
% fish swimming behavior, and food consumption. It uses semi-Lagrangian 
% advection to trace fish backward along the flow field and compute their 
% new location.
%
% INPUTS:
%   Predator            : [n x m] array of fish biomass  [g/m^3]
%   Prey                : [n x m] array of available food [g/m^3]
%   current             : [n x m x 2] with fields 'u' and 'v' for velocity 
%                         in x/y directions [m/s]
%                       : u = [ 1:m, 1:n, 1 ]
%                       : v = [ 1:m, 1:n, 2 ]
%   dt                  : timestep [s]
%   dx_m, dy_m          : grid spacing in x and y directions [m]
%   fish_speed          : scalar, max swimming speed of fish [m/s]
%   mask                : [n x m] logical array of valid ocean domain
%   m, n                : grid dimensions
%
% OUTPUTS:
%   Predator            : updated fish biomass field after advection 
%                         and feeding
%
% NOTES:
%   - The function handles boundary conditions using masking.
%   - The advection uses backward particle tracing to conserve fish biomass.
%
% AUTHOR: [JARED BRZENSKI
% DATE  : 30-06-2025
% -------------------------------------------------------------------------
%% The core of the advection fish code.
    % COBALT Inputs
    % mask off any advected food
    Prey = Prey .* mask;  
    %
    % Calculate the percent change between cells
    percent_more_food_full = percentMoreFoodFourWay( Prey );
    % Calculate the apparent current, taking swim speed and current
    apparent_current_full = ApparentCurrentFull( current, ...
                                                 fish_speed, ...
                                                 percent_more_food_full);
    %
    Fish_OG = Predator; % Hold original fish count
    %
    Flux = zeros(size(Predator));
    %
    for j=1:m
        for i=1:n
            idx = [i j];
            if Fish_OG(i,j) > 0
                
               Flux =  moveFishFunction( Predator, ...
                                         Flux, ...
                                         Prey, ...
                                         apparent_current_full, ...
                                         idx, ...
                                         dt, ...
                                         dx_m, ...
                                         dy_m );
            end
        end
    end
    Predator = Predator + Flux;
    %
end