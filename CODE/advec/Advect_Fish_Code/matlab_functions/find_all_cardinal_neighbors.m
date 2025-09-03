function neighbors_all = find_all_cardinal_neighbors(lat, lon)
% Find cardinal neighbors for all (i,j) in the grid
% Output: neighbors_all.direction(i,j,:) = [i_neigh, j_neigh]

    [ni, nj] = size(lat);

    % Preallocate
    directions = ["north", "south", "east", "west"];
    for d = directions
        neighbors_all.(d) = NaN(ni, nj, 2);
    end

    for j = 1:nj
        for i = 1:ni
            neighbors = find_cardinal_neighbors_tripolar(lat, lon, i, j);
            for d = directions
                neighbors_all.(d)(i,j,:) = neighbors.(d);
            end
        end
    end
end