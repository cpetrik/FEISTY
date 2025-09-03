function [i_north, j_north] = find_northern_neighbor(lat, lon, i, j)
% Find the physically closest northern neighbor of grid point (i,j)
% lat, lon: [ni x nj] matrices of latitudes and longitudes
% i, j: indices of the current grid cell
% Returns i_north, j_north: indices of the northern neighbor

    % Get the current point's coordinates
    lat0 = lat(i,j);
    lon0 = lon(i,j);

    % Flatten the lat/lon grid
    [ni, nj] = size(lat);
    latv = lat(:);
    lonv = lon(:);

    % Approximate angular distance (avoid needing Mapping Toolbox)
    % Prioritize only *northward* points
    north_mask = latv > lat0;

    % Use cosine-weighted Euclidean distance
    dist2 = (latv - lat0).^2 + (cosd(lat0).*(lonv - lon0)).^2;

    % Mask out non-northern points
    dist2(~north_mask) = inf;

    % Find the minimum distance among northern points
    [~, idx] = min(dist2);

    % Convert linear index to subscripts
    [i_north, j_north] = ind2sub([ni, nj], idx);
end