function [ Food ] = makeFood( varargin )
% Make food distribution
% scenario, 
% Food, 
% X
% Y
% x0, 
% y0, 
% amplitude, 
% concentration
%
numArgs = length(varargin);

switch numArgs
    case 4 % gradient in longer direction [0,1]
        scenario = varargin{1};
        Food = varargin{2};
        X = varargin{3};
        Y = varargin{4};

    case 8 % Gaussian Hump at (x0, y0)
        scenario = varargin{1};
        Food = varargin{2};
        X = varargin{3};
        Y = varargin{4};
        x0 = varargin{5};
        y0 = varargin{6};
        amplitude = varargin{7};
        concentration = varargin{8};

    otherwise
        error(' Unsupported number of args for makeFood');
end

[ n m ] = size( Food );

if strcmp( 'gradient', scenario )
    % Creating gradient (meshgrid) in X and Y
    if m > n % more columns, longer X direction
        scale = max(X(:));
        Food = 1 - X ./ scale .* ones(size(Food));
    else % more rows, longer Y direction
        scale = max(Y(:));
        Food = Y ./ scale .* ones(size(Food));
    end
    Food(:,1)   = 0.0;
    Food(:,end) = 0.0;
    Food(1,:)   = 0.0;
    Food(end,:) = 0.0;

elseif strcmp( 'slope', scenario )
    % Creating gradient (meshgrid) in X and Y
    if m > n % more columns, longer X direction
        scale = max(X(:));
        slope = 25;
        xmid = scale/2;
        Food = slope .* (X ./ scale) .* ones(size(Food));
    else % more rows, longer Y direction
        scale = max(Y(:));
        Food = Y ./ scale .* ones(size(Food));
    end
    Food(:,1)   = 0.0;
    Food(:,end) = 0.0;
    Food(1,:)   = 0.0;
    Food(end,:) = 0.0;

elseif strcmp( 'hump', scenario )
    % Gaussian distribution at (x0, y0), with a fixed concentration outside
    % the hump region
    sigma_x = 250; % Standard deviation in x
    sigma_y = 250; % Standard deviation in y
    initial_concentration = concentration; % average density
    % Define the Gaussian function
    Food = amplitude * exp(-((X-x0).^2 / (2*sigma_x^2) + (Y-y0).^2 / (2*sigma_y^2)));
    Food( Food < initial_concentration) = initial_concentration;

elseif strcmp( 'humpfish', scenario )
    % Gaussian distribution at (x0, y0), with a fixed concentration outside
    % the hump region
    sigma_x = 15; % Standard deviation in x
    sigma_y = 15; % Standard deviation in y
    initial_concentration = concentration; % average density
    % Define the Gaussian function
    Food = amplitude * exp(-((X-x0).^2 / (2*sigma_x^2) + (Y-y0).^2 / (2*sigma_y^2)));
    Food( Food < initial_concentration) = initial_concentration;

else
    error('Could not find scenario to use, returning Food');
end
    Food(:,1)   = 0.0;
    Food(:,end) = 0.0;
    Food(1,:)   = 0.0;
    Food(end,:) = 0.0;

end

        
