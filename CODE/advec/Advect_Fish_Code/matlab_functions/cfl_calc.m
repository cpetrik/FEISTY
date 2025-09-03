%% Grid lcoation test
%
%

close all; clear all; clc;

datafile='../data/feb152013_run25_ocean.198801-200712_uh200_vh200.mat';
gridfile='../data/Data_hindcast_grid_cp2D.mat';

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

%% CFL Calculations
dx_grid = GRD.dxtn;
dy_grid = GRD.dyte;

fast_numerator = 60*60*24*1;
slow_numerator = 60*60*24*0.1;

v_cfl = fast_numerator ./ dy_grid;
u_cfl = fast_numerator ./ dx_grid;

figure(11)
tiledlayout(2,3)
% Plot u_cfl
nexttile
v_cfl(v_cfl>2) = 2;
surf(geolat_t, geolon_t, v_cfl)
shading interp
colormap jet
hold on
mask = v_cfl;
mask(mask<=1) = NaN;
surf(geolat_t, geolon_t, mask, 'FaceColor', 'red', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('lat'); ylabel('lon'); colorbar;
caxis([0 2]);
title('V CFL for 1m/s max speed, dt = 1');
view(90, -90);

nexttile
v_cfl = fast_numerator ./ dy_grid ./ 2;
v_cfl(v_cfl>2) = 2;
surf(geolat_t, geolon_t, v_cfl)
shading interp
colormap jet
hold on
mask = v_cfl;
mask(mask<=1) = NaN;
surf(geolat_t, geolon_t, mask, 'FaceColor', 'red', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('lat'); ylabel('lon'); colorbar;
caxis([0 2]);
title('V CFL for 1m/s max speed, dt = 0.5 day');
view(90, -90);

nexttile
v_cfl = fast_numerator ./ dy_grid ./ 4;
v_cfl(v_cfl>2) = 2;
surf(geolat_t, geolon_t, v_cfl)
shading interp
colormap jet
hold on
mask = v_cfl;
mask(mask<=1) = NaN;
surf(geolat_t, geolon_t, mask, 'FaceColor', 'red', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('lat'); ylabel('lon'); colorbar;
caxis([0 2]);
title('V CFL for 1m/s max speed, dt = 0.25 day');
view(90, -90);

nexttile
u_cfl(u_cfl>2) = 2;
surf(geolat_t, geolon_t, u_cfl)
shading interp
colormap jet
hold on
mask = u_cfl;
mask(mask<=1) = NaN;
surf(geolat_t, geolon_t, mask, 'FaceColor', 'red', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('lat'); ylabel('lon'); colorbar;
caxis([0 2]);
title('U CFL for 1m/s max speed, dt = 1 day');
view(90,-90)

nexttile
u_cfl = fast_numerator ./ dx_grid ./ 2;

u_cfl(u_cfl>2) = 2;
surf(geolat_t, geolon_t, u_cfl)
shading interp
colormap jet
hold on
mask = u_cfl;
mask(mask<=1) = NaN;
surf(geolat_t, geolon_t, mask, 'FaceColor', 'red', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('lat'); ylabel('lon'); colorbar;
caxis([0 2]);
title('U CFL for 1m/s max speed, dt = 0.5 day');
view(90,-90)
v_cfl = slow_numerator ./ dy_grid;
u_cfl = slow_numerator ./ dx_grid;


nexttile
u_cfl = fast_numerator ./ dx_grid ./ 4;

u_cfl(u_cfl>2) = 2;
surf(geolat_t, geolon_t, u_cfl)
shading interp
colormap jet
hold on
mask = u_cfl;
mask(mask<=1) = NaN;
surf(geolat_t, geolon_t, mask, 'FaceColor', 'red', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('lat'); ylabel('lon'); colorbar;
caxis([0 2]);
title('U CFL for 1m/s max speed, dt = 0.25 day');
view(90,-90)
v_cfl = slow_numerator ./ dy_grid;
u_cfl = slow_numerator ./ dx_grid;


%% Slow speed
figure(13)
tiledlayout(2,2)
nexttile
v_cfl(v_cfl>2) = 2;
surf(geolat_t, geolon_t, v_cfl)
shading interp
colormap jet
hold on
mask = v_cfl;
mask(mask<=1) = NaN;
surf(geolat_t, geolon_t, mask, 'FaceColor', 'red', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('lat'); ylabel('lon'); colorbar;
caxis([0 2]);
title('V CFL for 0.1m/s max speed dt = 1 day');
view(90, -90);

nexttile
v_cfl = slow_numerator ./ dy_grid ./ 2;
v_cfl(v_cfl>2) = 2;
surf(geolat_t, geolon_t, v_cfl)
shading interp
colormap jet
hold on
mask = v_cfl;
mask(mask<=1) = NaN;
surf(geolat_t, geolon_t, mask, 'FaceColor', 'red', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('lat'); ylabel('lon'); colorbar;
caxis([0 2]);
title('V CFL for 0.1m/s max speed dt = 0.5 days');
view(90, -90);

nexttile
u_cfl(u_cfl>2) = 2;
surf(geolat_t, geolon_t, u_cfl)
shading interp
colormap jet
hold on
mask = u_cfl;
mask(mask<=1) = NaN;
surf(geolat_t, geolon_t, mask, 'FaceColor', 'red', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('lat'); ylabel('lon'); colorbar;
caxis([0 2]);
title('U CFL for 0.1m/s max speed dt = 1 day');
view(90,-90)

nexttile
u_cfl = slow_numerator ./ dx_grid ./ 2;
u_cfl(u_cfl>2) = 2;
surf(geolat_t, geolon_t, u_cfl)
shading interp
colormap jet
hold on
mask = u_cfl;
mask(mask<=1) = NaN;
surf(geolat_t, geolon_t, mask, 'FaceColor', 'red', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('lat'); ylabel('lon'); colorbar;
caxis([0 2]);
title('U CFL for 0.1m/s max speed dt = 0.5 days');
view(90,-90)