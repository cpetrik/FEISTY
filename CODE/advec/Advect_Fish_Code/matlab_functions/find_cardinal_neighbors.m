function neighbors = find_cardinal_neighbors(lat, lon, i, j)
% Find the physically closest N, S, E, W neighbors of grid point (i,j)
% lat, lon: [ni x nj] latitude and longitude matrices
% Returns neighbors.(north/south/east/west) = [i, j] or [NaN, NaN]

    [ni, nj] = size(lat);
    lat0 = lat(i,j);
    lon0 = lon(i,j);

    % Flatten lat/lon grid
    latv = lat(:);
    lonv = lon(:);

    % Exclude self
    latv(i + (j-1)*ni) = NaN;
    lonv(i + (j-1)*ni) = NaN;

    % Compute bearing (azimuth) from (lat0, lon0) to each point
    az = compute_azimuth(lat0, lon0, latv, lonv);  % degrees from north

    % Approximate distance
    dist2 = (latv - lat0).^2 + (cosd(lat0).*(lonv - lon0)).^2;

    % Direction masks (±22.5° around cardinal)
    dir_masks.north = (az >= 337.5 | az < 22.5);
    dir_masks.east  = (az >= 67.5  & az < 112.5);
    dir_masks.south = (az >= 157.5 & az < 202.5);
    dir_masks.west  = (az >= 247.5 & az < 292.5);

    % Find closest neighbor in each direction
    neighbors = struct();
    for dir = ["north", "south", "east", "west"]
        mask = dir_masks.(dir);
        d = dist2;
        d(~mask) = inf;
        [~, idx] = min(d);
        if isinf(d(idx))
            neighbors.(dir) = [NaN, NaN];
        else
            [ii, jj] = ind2sub([ni, nj], idx);
            neighbors.(dir) = [ii, jj];
        end
    end
end

function az = compute_azimuth(lat1, lon1, lat2, lon2)
% Computes azimuth from (lat1, lon1) to (lat2, lon2) in degrees
% Using spherical law of cosines formula
    dlon = deg2rad(lon2 - lon1);
    lat1 = deg2rad(lat1);
    lat2 = deg2rad(lat2);

    x = cos(lat2) .* sin(dlon);
    y = cos(lat1) .* sin(lat2) - sin(lat1) .* cos(lat2) .* cos(dlon);

    az = atan2d(x, y);
    az = mod(az, 360);  % ensure range [0, 360)
end
