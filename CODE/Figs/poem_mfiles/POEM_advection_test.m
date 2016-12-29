clear all
nc64startup

Uth_200 = nc_varget('/Volumes/GFDL/GCM_DATA/CORE-forced/feb152013_run25_ocean.198801-200712_uh200_vh200.nc','Uth_200',[0 0 0],[1 200 360]);
Vth_200 = nc_varget('/Volumes/GFDL/GCM_DATA/CORE-forced/feb152013_run25_ocean.198801-200712_uh200_vh200.nc','Vth_200',[0 0 0],[1 200 360]);

geolon_t = nc_varget('/Volumes/GFDL/GCM_DATA/Hindcast/grid_spec.nc','geolon_t',[0 0],[200 360]);
geolat_t = nc_varget('/Volumes/GFDL/GCM_DATA/Hindcast/grid_spec.nc','geolat_t',[0 0],[200 360]);
dxtn = nc_varget('/Volumes/GFDL/GCM_DATA/Hindcast/grid_spec.nc','dxtn',[0 0],[200 360]);
dyte = nc_varget('/Volumes/GFDL/GCM_DATA/Hindcast/grid_spec.nc','dyte',[0 0],[200 360]);
ht = nc_varget('/Volumes/GFDL/GCM_DATA/Hindcast/grid_spec.nc','ht',[0 0],[200 360]);
area = nc_varget('/Volumes/GFDL/GCM_DATA/Hindcast/grid_spec.nc','AREA_OCN',[0 0],[200 360]);

% rotate everything so that the first dimension is longitudes, w/1
% corresponding the the western-most point on the grid and moving from west
% to east; the second dimension is latitude with 1 corresponding to
% Antarctica and moving north;

Uth_200 = flipud(rot90(Uth_200));
Vth_200 = flipud(rot90(Vth_200));
geolon_t = flipud(rot90(geolon_t));
geolat_t = flipud(rot90(geolat_t));
dxtn = flipud(rot90(dxtn));
dyte = flipud(rot90(dyte));
ht = flipud(rot90(ht));
area = flipud(rot90(area))*510072000*1e6;
area = max(area,1);

% depth of the surface layer, 200m or less
eps = 1;
dep = min(ht,200);
dep = max(dep,eps);

%define a patch to advect
TF = zeros(360,200);
TF2 = zeros(360,200);
lonmin = -280;
lonmax = 80;
latmin = -90;
latmax = 90;
aa = find( (geolon_t > lonmin) & (geolon_t < lonmax) & (geolat_t > latmin) & ...
           (geolat_t < latmax) & (ht > 0) );
TF(aa) = 1;
total_mass(1) = sum(TF(:).*area(:));

% Following Advect_upwind_2D
dt = 3600;
ni = 360;
nj = 200;
isd = 1;
jsd = 1;
ied = ni;
jed = nj;

ntime = 365*24;
uvel = Uth_200;
vvel = Vth_200;

mask = zeros(ni,nj);
aa = find(ht > 0);
mask(aa) = 1;

fe = zeros(ni,nj);
fn = zeros(ni,nj);

% Advection loop
for n = 1:ntime
n
% Westward flux
for j = jsd:jed
    for i = isd:ied
        velocity = 0.5*uvel(i,j);
        upos = velocity + abs(velocity);
        uneg = velocity - abs(velocity);

        % define only for ocean cells
        if (mask(i,j) > 0)
        
        if (i == ied)
            fe(i,j) = dyte(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(isd,j))/dep(isd,j)* ...
                mask(i,j)*mask(isd,j);
        else
            fe(i,j) = dyte(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(i+1,j)/dep(i+1,j))* ...
                mask(i,j)*mask(i+1,j);
        end
        
        end
    end
end

% northward flux
for j = jsd:jed
    for i = isd:ied
        velocity = 0.5*vvel(i,j);
        upos = velocity + abs(velocity);
        uneg = velocity - abs(velocity);
            
        if (j < jed)
            fn(i,j) = dxtn(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(i,j+1)/dep(i,j+1))* ...
                mask(i,j)*mask(i,j+1);
        else
            fn(i,j) = dxtn(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(ni-i+1,j)/dep(ni-i+1,j))* ...
                mask(i,j)*mask(ni-i+1,j);
        end
        
    end
end

% combine fluxes
for j = jsd:jed
    for i = isd:ied
        if (j > 1)
            if (i > 1)
                upwind(i,j) = mask(i,j).*(fe(i-1,j)-fe(i,j)+fn(i,j-1)-fn(i,j));
            else
                upwind(i,j) = mask(i,j).*(fe(ied,j)-fe(i,j)+fn(i,j-1)-fn(i,j));
            end
        end
    end
end

% update tracers
for j = jsd:jed
    for i = isd:ied
        TF2(i,j) = TF(i,j) + (dt*upwind(i,j))/area(i,j);
    end
end
total_mass(n+1) = sum(TF2(:).*area(:));

% plot, do mass balance, reset tracer fields
aa = find(ht == 0);
TF(aa) = -999;
if n == 1
    figure(1)
    surf(geolon_t,geolat_t,TF); view(2); shading interp; caxis([0 1]);
    %pause
end

TF2(aa) = -999;

if n == 8760
figure(1)
clf
surf(geolon_t,geolat_t,TF2); view(2); shading interp; caxis([0 1]);
pdiff = 100*(total_mass(n+1) - total_mass(n))/total_mass(n);
title(['%diff = ', num2str(pdiff,'%10.8f')]);
pause
end

TF = TF2;

end



        








