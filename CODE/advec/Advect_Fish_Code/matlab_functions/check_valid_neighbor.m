function valid = check_valid_neighbor(neighbors, mask)
    % Initialize output
    valid = false(4,1);

    % Get mask size
    [nRows, nCols] = size(mask);

    for k = 1:4
        i = neighbors(k, 1);
        j = neighbors(k, 2);

        % Check if indices are within bounds
        if i >= 1 && i <= nRows && j >= 1 && j <= nCols
            % Check if mask value is valid (i.e., 1)
            if mask(i, j) == 1
                valid(k) = true;
            end
        end
    end
end