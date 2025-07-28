function neighbors = get_neighbors_from_struct(neighbors_all, i, j)
% Extracts the [i, j] neighbor indices for a specific point (i,j)
% from the full neighbors_all struct.

    directions = ["north", "south", "east", "west"];
    for d = directions
        neighbors.(d) = squeeze(neighbors_all.(d)(i,j,:))';
    end
end