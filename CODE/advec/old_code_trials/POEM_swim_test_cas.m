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

%%
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

%% Swimming behavior
wgt = 2.5; 	%M=2.5; L=2500.0
T = 15.0;
w = ((3.9*wgt.^0.13 * exp(0.149*T)) /100);
Q = zeros(360,200);
Q(aa) = w;

%define value to maximize
nu = ht; %go to deep
%nu = -1.0 * ht; %go to shallow

%% Following Advect_upwind_2D
ni = 360;
nj = 200;
isd = 1;
jsd = 2;
ied = ni;
jed = nj;

dt = 3600;
ntime = 365*24;

mask = zeros(ni,nj);
aa = find(ht > 0);
mask(aa) = 1;

fe = zeros(ni,nj);
fn = zeros(ni,nj);

%% Advection loop
for n = 1:ntime
    n
    % Find desired cell
    KK = zeros(360,200);
    for j=jsd:jed
        for i=isd:ied
            if (j==nj)
                if (i==1)
                    [temp,KK(i,j)] = max([nu(i,j),nu(ni-i+1,j),nu(i,j-1),nu(ied,j),nu(i+1,j)]);
                elseif (i==ni)
                    [temp,KK(i,j)] = max([nu(i,j),nu(ni-i+1,j),nu(i,j-1),nu(i-1,j),nu(isd,j)]);
                else
                    [temp,KK(i,j)] = max([nu(i,j),nu(ni-i+1,j),nu(i,j-1),nu(i-1,j),nu(i+1,j)]);
                end
            else
                if (i==1)
                    [temp,KK(i,j)] = max([nu(i,j),nu(i,j+1),nu(i,j-1),nu(ied,j),nu(i+1,j)]);
                elseif (i==ni)
                    [temp,KK(i,j)] = max([nu(i,j),nu(i,j+1),nu(i,j-1),nu(i-1,j),nu(isd,j)]);
                else
                    [temp,KK(i,j)] = max([nu(i,j),nu(i,j+1),nu(i,j-1),nu(i-1,j),nu(i+1,j)]);
                end
            end
        end %i
    end %j
    KK = KK .* mask;
    
    % Adjust velocity towards cell
    Ua = zeros(360,200);
    Va = zeros(360,200);
    I1 = find(KK == 1);
    I2 = find(KK == 2);
    I3 = find(KK == 3);
    I4 = find(KK == 4);
    I5 = find(KK == 5);
    Va(I2) = Va(I2) + Q(I2);
    Va(I3) = Va(I3) - Q(I3);
    Ua(I4) = Ua(I4) - Q(I4);
    Ua(I5) = Ua(I5) + Q(I5);
    
    % Westward flux
    for j = jsd:jed
        for i = isd:ied
            velocity = 0.5*Ua(i,j);
            upos = velocity + abs(velocity);
            uneg = velocity - abs(velocity);
            
            % define only for ocean cells
            if (mask(i,j) > 0)
                
                if (i == ied)
                    %fe(i,j) = dyte(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(isd,j)/dep(isd,j))* ...
                    %mask(i,j)*mask(isd,j);
                    fe(i,j) = dyte(i,j)*(upos.*TF(i,j) + uneg.*TF(isd,j))* ...
                        mask(i,j)*mask(isd,j);
                else
                    %fe(i,j) = dyte(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(i+1,j)/dep(i+1,j))* ...
                    %mask(i,j)*mask(i+1,j);
                    fe(i,j) = dyte(i,j)*(upos.*TF(i,j) + uneg.*TF(i+1,j))* ...
                        mask(i,j)*mask(i+1,j);
                end
                
            end
        end
    end
    
    % northward flux
    for j = jsd:jed
        for i = isd:ied
            velocity = 0.5*Va(i,j);
            upos = velocity + abs(velocity);
            uneg = velocity - abs(velocity);
            
            if (j < jed)
                %fn(i,j) = dxtn(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(i,j+1)/dep(i,j+1))* ...
                %mask(i,j)*mask(i,j+1);
                fn(i,j) = dxtn(i,j)*(upos.*TF(i,j) + uneg.*TF(i,j+1))* ...
                    mask(i,j)*mask(i,j+1);
            else
                %fn(i,j) = dxtn(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(ni-i+1,j)/dep(ni-i+1,j))* ...
                %mask(i,j)*mask(ni-i+1,j);
                fn(i,j) = dxtn(i,j)*(upos.*TF(i,j) + uneg.*TF(ni-i+1,j))* ...
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
    
    TF = TF2;
    
end

figure(2)
clf
surf(geolon_t,geolat_t,TF2); view(2); shading interp; caxis([0 1]);
pdiff = 100*(total_mass(n+1) - total_mass(n))/total_mass(n);
title(['%diff = ', num2str(pdiff,'%10.3e')]);
%pause

figure
surf(geolon_t,geolat_t,fe);
view(2); shading interp;
colorbar;
title('fe')

figure
surf(geolon_t,geolat_t,fn);
view(2); shading interp;
colorbar;
title('fn')