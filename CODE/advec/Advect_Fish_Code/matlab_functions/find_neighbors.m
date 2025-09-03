%% Grid lcoation test
%
%

close all; clear all; clc;

datafile='../data/feb152013_run25_ocean.198801-200712_uh200_vh200.mat';
gridfile='../data/Data_hindcast_grid_cp2D.mat';

addpath(genpath('neighbors_functions'))

load(datafile);
load(gridfile);

% Assume lat and lon are [ni x nj]
[nx, ny] = size(geolat_t);

nbs = nan(nx, ny, 4);  % 4 neighbors: N, S, E, W

for i = 1:nx
    for j = 1:ny
        if GRD.mask(i,j)  % use wet(i,j) or mask(i,j) if available

            % North
            if i > 1 && GRD.mask(i-1,j)
                nbs(i,j,1) = sub2ind([nx,ny], i-1, j);
            end

            % South
            if i < nx && GRD.mask(i+1,j)
                nbs(i,j,2) = sub2ind([nx,ny], i+1, j);
            end

            % East (wraparound)
            jj = mod(j, ny) + 1;
            if GRD.mask(i,jj)
                nbs(i,j,3) = sub2ind([nx,ny], i, jj);
            end

            % West (wraparound)
            jj = mod(j-2, ny) + 1;
            if GRD.mask(i,jj)
                nbs(i,j,4) = sub2ind([nx,ny], i, jj);
            end
        end
    end
end


x_vec = linspace(1,nx,nx);
y_vec = linspace(1,ny,ny);
[XX,YY] = meshgrid(y_vec, x_vec);
xpoint = 120;
ypoint = 200;

nbs(xpoint,ypoint,:)

temp = zeros(nx,ny);

%geolon_t(geolon_t < -180) = geolon_t(geolon_t < -180) + 360;

[i_north, j_north] = find_northern_neighbor(geolat_t, geolon_t, xpoint, ypoint);
neighbors = find_cardinal_neighbors( geolat_t, geolon_t, xpoint, ypoint);


tic
% Time consuming process, but faster to do it now than individually in a
% loop. Takes about 2 minutes with 72000 points
all_neighbors = find_all_cardinal_neighbors( geolat_t, geolon_t );
toc

figure(1)
scatter(geolon_t, geolat_t); hold on;
plot(geolon_t(xpoint, ypoint), geolat_t(xpoint, ypoint), 'ko','LineWidth',5);
plot(geolon_t(neighbors.east(1), neighbors.east(2)), geolat_t(neighbors.east(1), neighbors.east(2)), 'mo','LineWidth',3);
plot(geolon_t(neighbors.west(1), neighbors.west(2)), geolat_t(neighbors.west(1), neighbors.west(2)), 'go','LineWidth',3);
plot(geolon_t(neighbors.north(1), neighbors.north(2)), geolat_t(neighbors.north(1), neighbors.north(2)), 'ro','LineWidth',3);
plot(geolon_t(neighbors.south(1), neighbors.south(2)), geolat_t(neighbors.south(1), neighbors.south(2)), 'bo','LineWidth',3);


alt_n = get_neighbors_from_struct( all_neighbors, xpoint, ypoint );

figure(2)
scatter(geolon_t, geolat_t); hold on;
plot(geolon_t(xpoint, ypoint), geolat_t(xpoint, ypoint), 'ko','LineWidth',5);
plot(geolon_t(alt_n.east(1), alt_n.east(2)), geolat_t(alt_n.east(1), alt_n.east(2)), 'mo','LineWidth',3);
plot(geolon_t(alt_n.west(1), alt_n.west(2)), geolat_t(alt_n.west(1), alt_n.west(2)), 'go','LineWidth',3);
plot(geolon_t(alt_n.north(1), alt_n.north(2)), geolat_t(alt_n.north(1), alt_n.north(2)), 'ro','LineWidth',3);
plot(geolon_t(alt_n.south(1), alt_n.south(2)), geolat_t(alt_n.south(1), alt_n.south(2)), 'bo','LineWidth',3);

xpoint=100; ypoint = 100;
alt_n = get_neighbors_from_struct( all_neighbors, xpoint, ypoint );

figure(3)
scatter(geolon_t, geolat_t); hold on;
plot(geolon_t(xpoint, ypoint), geolat_t(xpoint, ypoint), 'ko','LineWidth',5);
plot(geolon_t(alt_n.east(1), alt_n.east(2)), geolat_t(alt_n.east(1), alt_n.east(2)), 'mo','LineWidth',3);
plot(geolon_t(alt_n.west(1), alt_n.west(2)), geolat_t(alt_n.west(1), alt_n.west(2)), 'go','LineWidth',3);
plot(geolon_t(alt_n.north(1), alt_n.north(2)), geolat_t(alt_n.north(1), alt_n.north(2)), 'ro','LineWidth',3);
plot(geolon_t(alt_n.south(1), alt_n.south(2)), geolat_t(alt_n.south(1), alt_n.south(2)), 'bo','LineWidth',3);
