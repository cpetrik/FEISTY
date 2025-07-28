function [ dir ] = MoreFoodFourWay( Food, idx, neighbors )
%% This finds food in four directions, and only looks for MORE food!
%
% INPUTS:
%   Food  : [m x n x 4] array of percent more/less food concentrations
%   idx   : index [m_i, n_i] of the specific location
%
% OUTPUT:
%   dir   : [4 x 1] logical array of where there is more food
%           in the four cardinal directions, up, down, left, and right
%
%
% Author: JARED BRZENSKI
% Date  : 30-06-2025
%
%  dir = MORE food than current cell = [ U D R L ]
% -------------------------------------------------------------------------
%% Core percent more food code
% Finds the direction of more food ONLY! 
    dir = [ 0 0 0 0 ];
    
    percent_more_food = [ 0 0 0 0 ];

    i = idx(1);
    j = idx(2);

    current_val = Food(i,j);

    [n, m] = size(Food);
    
    % Define neighbors relative to (i, j)
    % neighbors = [-1, 0;  % Up
    %               1, 0;  % Down
    %               0, 1;  % Right
    %               0, -1];% Left

    % Loop through each direction ( find more food )
    for k = 1:4
        %ni = i + neighbors(k, 1); % Neighbor row
        %nj = j + neighbors(k, 2); % Neighbor column
        
        ni = neighbors(k, 1);
        nj = neighbors(k, 2);

        % Check bounds and compare food values
        %if ni >= 1 && ni <= n && nj >= 1 && nj <= m % Valid neighbor
        if Food(ni, nj) > current_val
            dir(k) = 1;
        end
        %end
    end
    
end
