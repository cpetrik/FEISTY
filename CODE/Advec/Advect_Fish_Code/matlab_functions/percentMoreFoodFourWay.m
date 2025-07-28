function [ percent_more_food_full ] = percentMoreFoodFourWay( Food, neighbors_all )
%% This finds food in four directions, and then calculates percentage
%  more food in that direction.
%
% INPUTS:
%   Food    : [m x n] array of prey concentrations [unitless]
%
% OUTPUT:
%   precent_more_food_full : [ m x n x 4] percent more food for each cell
%                                         in the four cardinal directions, 
%                                         up, down, left, and right
%
% TEST : test with the function test_percentMoreFoodFourWay is avail.
%
% Author: JARED BRZENSKI
% Date  : 30-06-2025
%
%  dir = percent MORE than current cell = [ U D R L ]
% -------------------------------------------------------------------------
%% Core percent more food code
    % Define neighbors relative to (i, j)
    neighbors = [-1, 0;  % Up
                  1, 0;  % Down
                  0, 1;  % Right
                  0, -1];% Left

    [n, m] = size(Food);

    percent_more_food_full = zeros(n,m,4);

    percent_more_food = [ 0 0 0 0 ];

    for j=1:m
        for i=1:n
            current_val = Food(i,j);
            % Reset food vector each loop
            percent_more_food = [ 0 0 0 0 ];
            % Local Neighborhood
            neighborhood = get_neighbors_from_struct(neighbors_all, i, j);
            neighbors = [neighborhood.north; ...
                         neighborhood.south; ...
                         neighborhood.east; ...
                         neighborhood.west];
            % Loop through each direction ( find more food )
            for k = 1:4
                % ni = i + neighbors(k, 1); % Neighbor row
                % nj = j + neighbors(k, 2); % Neighbor column
                ni = neighbors(k, 1);
                nj = neighbors(k, 2);
                
                % Check bounds and compare food values
                %if ni >= 1 && ni <= n && nj >= 1 && nj <= m % Valid neighbor
                if Food(ni, nj) > current_val
                    percent_more_food(k) = min( Food(ni,nj)/current_val, 2 ) - 1;
                end
                %end
            end
            percent_more_food_full(i,j,:) = percent_more_food;
        end
    end
end
