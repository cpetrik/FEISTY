function neighbors = find_cardinal_neighbors_tripolar(lat, lon, i, j)
    [nRows, nCols] = size(lat);
    lat0 = lat(i,j);
    lon0 = lon(i,j);

    % Flatten the grid
    latVec = lat(:);
    lonVec = lon(:);

    current_flat_idx = sub2ind(size(lat), i, j);

    % Compute great-circle bearings from (lat0, lon0) to all other points
    bearings = compute_bearing(lat0, lon0, latVec, lonVec);  % in degrees
    dists = distance_gc(lat0, lon0, latVec, lonVec);         % in degrees or km

    % Define cardinal direction targets
    targets = struct('north', 0, 'east', 90, 'south', 180, 'west', -90);
    directions = fieldnames(targets);

    neighbors = struct();
    for k = 1:numel(directions)
        dir = directions{k};
        target_bearing = targets.(dir);

        % Find closest bearing match in direction
        angle_diff = angular_difference(bearings, target_bearing);
        angle_diff(current_flat_idx) = Inf;         % Exclude self
        [~, idx] = min(angle_diff + 0.01 * dists);  % prefer closer if angle is similar

        % Convert flat index back to (i,j)
        [ii, jj] = ind2sub(size(lat), idx);
        neighbors.(dir) = [ii, jj];
        %neighbors.(dir).index = [ii, jj];
        %neighbors.(dir).lat = lat(ii,jj);
        %neighbors.(dir).lon = lon(ii,jj);
    end
end


function b = compute_bearing(lat1, lon1, lat2, lon2)
    % All inputs in degrees
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);

    dLon = lon2 - lon1;
    x = sin(dLon) .* cos(lat2);
    y = cos(lat1) .* sin(lat2) - sin(lat1) .* cos(lat2) .* cos(dLon);
    b = atan2d(x, y);  % bearing in degrees
end

function d = distance_gc(lat1, lon1, lat2, lon2)
    % Great-circle distance using haversine formula
    R = 6371;  % km
    dlat = deg2rad(lat2 - lat1);
    dlon = deg2rad(lon2 - lon1);
    lat1 = deg2rad(lat1);
    lat2 = deg2rad(lat2);

    a = sin(dlat/2).^2 + cos(lat1).*cos(lat2).*sin(dlon/2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    d = R * c;
end

function diff = angular_difference(a1, a2)
    % Smallest absolute angular difference (in degrees)
    diff = mod(a1 - a2 + 180, 360) - 180;
    diff = abs(diff);
end
