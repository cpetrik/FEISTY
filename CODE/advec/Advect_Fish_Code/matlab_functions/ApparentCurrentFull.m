function [ dir_swim_speed_full ] = ApparentCurrentFull( current, fish_speed, percent_more_food )
%% Calculates apparent current due to fish possibly swimming with 
% or against the current. Units should be the same of fish_speed
% and current
%
% INPUTS:
%  
%   current             : [m x n x 2] with fields 'u' and 'v' for velocity
%                         in x/y directions [m/s]
%                       : u = [ 1:m, 1:n, 1 ]
%                       : v = [ 1:m, 1:n, 2 ]
%  fish_speed           : scalar, max swimming speed of fish [m/s]
%  percent_more_food    : [m x n x 4] array of food percentages in the four
%                         cardinal direction, up, down, right, left
%
% OUTPUT:
%   dir_swim_speed_full : [m x n x 4] array of swimming speeds, corrected
%                         for current, in the four cardinal direction 
%                         up, down, right, left.
%
% NOTES: 
%       - This array takes percent more food, in case we want to use a 
%         gradient for the swimming speed. We can adjust that here.
%
% AUTHOR: JARED BRZENSKI
% DATE  : 30-06-2025
% -------------------------------------------------------------------------
%% Core Apparent Current Full Code
    % Directions matrix
    directions = [
                -1, 0;  % Up
                 1, 0;  % Down
                 0, 1;  % Right
                 0, -1];% Left

    [n, m, o] = size(current);

    dir_swim_speed_full = zeros(n,m,4);

    for i=1:n
        for j=1:m
   
            U = current(i,j,1);
            V = current(i,j,2);
            
            dir_swim_speed_single = zeros(4,1);
            
            % Advection velocity, regardless of food
            
            % How fast to swim in each direction, based on percent more food

            per_more_food = reshape(percent_more_food(i,j,:), 4,1);

            % fourWayFishSpeed = mapSpeed( 'linear', per_more_food, fish_speed);
            unit_vector = zeros(size(per_more_food));
            unit_vector( per_more_food > 0 ) = 1;
            fourWayFishSpeed = fish_speed * unit_vector;
            
            % Total speed, current PLUS swimming speed
            % dir_swim_speed = vel + fourWayFishSpeed;
            dir_swim_speed_single = [ max( fourWayFishSpeed(1) + V, 0 ) ...
                               max( fourWayFishSpeed(2) - V, 0 ) ...
                               max( fourWayFishSpeed(3) + U, 0 ) ...
                               max( fourWayFishSpeed(4) - U, 0 ) ];

            dir_swim_speed_full(i,j,:) = dir_swim_speed_single;

        end
    end
end