clear all
close all

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'AREA_OCN','dat','dxtn','dyte','ht','geolon_t','geolat_t');
area = AREA_OCN;
area = (area)*510072000*1e6;
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

ni = 360;
nj = 200;
isd = 1;
jsd = 2;
ied = ni;
jed = nj;

%% Swimming behavior

load('/Volumes/GFDL/GCM_DATA/CORE-forced/ocean_cobalt.feb15_run25.1988-2007_phyt_zoop.mat',...
    'nmdz_avg200_88','nlgz_avg200_88')
load('/Volumes/GFDL/GCM_DATA/CORE-forced/ocean.186101-200512.temp_100_avg.mat',...
    'TEMP_100')

zoop = nmdz_avg200_88(:,:,end-11:end) + nlgz_avg200_88(:,:,end-11:end);
T100 = TEMP_100(:,:,end-11:end) - 273.15;

mid = find(T100 == min(T100(:)));
T100(mid) = NaN;
zid = find(zoop == min(zoop(:)));
zoop(zid) = NaN;

jdmo = 15:30:365;
zoopi = zeros(ni,nj,365);
ti = zoopi;
for i=1:ni
    for j=1:nj
        zoopi(i,j,:) = interp1(jdmo,squeeze(zoop(i,j,:)),1:365,'linear','extrap');
        ti(i,j,:) = interp1(jdmo,squeeze(T100(i,j,:)),1:365,'linear','extrap');
    end
end

clear TEMP_100 nmdz_avg200_88 nlgz_avg200_88

% L_m=10^((log10(20)+log10(200))/2);
% w = exp(0.063*(T100-15.0)) * 0.5*L_m*1e-3;
% Q = zeros(360,200);
% Q(aa) = w(aa);

wgt = 2.5;
T = 15.0;
w = ((3.9*wgt.^0.13 * exp(0.149*T)) /100);
Q = zeros(360,200);
Q(aa) = w;

%define value to maximize
%nu = ht; %go to deep
%nu = -1.0 * ht; %go to shallow
nu = zoop;

%% Following Advect_upwind_2D

dt = 3600;
ntime = 365*24;

mask = zeros(ni,nj);
aa = find(ht > 0);
mask(aa) = 1;

fe = zeros(ni,nj);
fn = zeros(ni,nj);

%% Advection loop
n=0;
%for n = 1:ntime
for d=1:365
    nu = zoopi(:,:,d);
    for h=1:24
        n=n+1
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
                        %                     fe(i,j) = dyte(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(isd,j)/dep(isd,j))* ...
                        %                         mask(i,j)*mask(isd,j);
                        fe(i,j) = dyte(i,j)*(upos.*TF(i,j) + uneg.*TF(isd,j))* ...
                            mask(i,j)*mask(isd,j);
                    else
                        %                     fe(i,j) = dyte(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(i+1,j)/dep(i+1,j))* ...
                        %                         mask(i,j)*mask(i+1,j);
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
                    %                 fn(i,j) = dxtn(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(i,j+1)/dep(i,j+1))* ...
                    %                     mask(i,j)*mask(i,j+1);
                    fn(i,j) = dxtn(i,j)*(upos.*TF(i,j) + uneg.*TF(i,j+1))* ...
                        mask(i,j)*mask(i,j+1);
                else
                    %                 fn(i,j) = dxtn(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(ni-i+1,j)/dep(ni-i+1,j))* ...
                    %                     mask(i,j)*mask(ni-i+1,j);
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
end
%
figure(2)
clf
surf(geolon_t,geolat_t,TF2); view(2); shading interp; caxis([0 1]);
pdiff = 100*(total_mass(n+1) - total_mass(n))/total_mass(n);
title(['%diff = ', num2str(pdiff,'%10.3e')]);
%pause

figure(3)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,TF2);
shading interp; caxis([0 1]);
colorbar
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');

figure(4)
plot(log(total_mass),'linewidth',2)

tdiff = diff(total_mass) ./ total_mass(1:end-1);
figure(5)
plot(100*tdiff,'linewidth',2)

% for t=1:75:365
%     figure
%     surf(geolon_t,geolat_t,zoopi(:,:,t)); view(2); shading interp; 
%     %caxis([0 1]);
%     title(['zoop ' num2str(t)]);
% end

save('/Volumes/GFDL/CSV/advect_tests/matbio_2Dadvect_swim_zoop2_test_global_vel0_dt1hr')

%%
eps = 1e-9;

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

%%
figure
m_proj('stereographic','lat',90,'long',30,'radius',20);
m_pcolor(geolon_t,geolat_t,fe);
shading interp
colorbar
colormap('jet')
title('fe')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');

figure
m_proj('stereographic','lat',90,'long',30,'radius',20);
m_pcolor(geolon_t,geolat_t,fn);
shading interp
colorbar
colormap('jet')
title('fn')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');

%%
figure
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,dxtn);
shading interp
colorbar
colormap('jet')
title('dxtn')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');

figure
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,dyte);
shading interp
colorbar
colormap('jet')
title('dyte')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');


%%
figure
bar(fn(:,200))
id = find(fn(:,200)<-2e5)

test=zeros(360,200);
test(id,200)=1;

figure
m_proj('stereographic','lat',90,'long',30,'radius',20);
m_pcolor(geolon_t,geolat_t,test);
shading interp
title('weird grid cells')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');

%%
figure
m_proj('stereographic','lat',90,'long',30,'radius',20);
m_pcolor(geolon_t,geolat_t,KK);
shading flat
title('KK')
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');








