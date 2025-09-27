function Zsmooth = smooth2nan_conservative_mask(Z, window)
    % Truly conservative, NaN-preserving smoothing
    % Redistributes mass only within the valid region.
    %
    % Z: 2D array with NaNs
    % window: size of square window (e.g., 3 for 3x3)

    if nargin < 2
        window = 3;
    end

    mask = ~isnan(Z);    % valid data mask, basically where nans are.
    Zfilled = Z;
    Zfilled(~mask) = 0;

    kernel = ones(window);

    % Compute sum of neighbors and count of valid neighbors
    sumValues   = conv2(Zfilled, kernel, 'same');
    countValues = conv2(mask,   kernel, 'same');

    % Compute normalized smoothing weights for each cell
    % Instead of a simple average, normalize so that
    % the sum of contributions from neighbors = original value.
    Zsmooth = zeros(size(Z));
    for j = 1:size(Z,2)
        for i = 1:size(Z,1)
            if mask(i,j)
                % Get local window indices
                i1 = max(1, i-floor(window/2));
                i2 = min(size(Z,1), i+floor(window/2));
                j1 = max(1, j-floor(window/2));
                j2 = min(size(Z,2), j+floor(window/2));

                localMask = mask(i1:i2, j1:j2);
                localValues = Zfilled(i1:i2, j1:j2);

                % Compute weights only over valid neighbors
                w = localMask;
                w = w / sum(w(:));  % normalize to 1

                % Distribute mass from (i,j) proportionally
                Zsmooth(i1:i2, j1:j2) = Zsmooth(i1:i2, j1:j2) + ...
                                         Z(i,j) * w;
            end
        end
    end

    % Apply mask to keep NaNs
    Zsmooth(~mask) = NaN;
end
