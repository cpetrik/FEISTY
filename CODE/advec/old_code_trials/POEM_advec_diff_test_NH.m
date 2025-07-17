% Advection and diffusion using upwind scheme adapted from MOM code

clear all
close all

fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/advect_tests/';

load('/Volumes/GFDL/GCM_DATA/CORE-forced/Natasha_jellies/ocean.194801-200712_uh50_vh50.mat',...
    'uh50','vh50','u50','v50');

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'AREA_OCN','dat','dxtn','dyte','ht','geolon_t','geolat_t');

cname='AtlNH_dt1hr_velH50_Jan88_b100_v2';

%%
area = AREA_OCN;
area = (area)*510072000*1e6;
area = max(area,1);

% grid size 
[ni,nj] = size(geolon_t);
isd = 1;
jsd = 1;
ied = ni;
jed = nj;

% depth of the surface layer, 50m or less
eps = 1;
dep = min(ht,50);
dep = max(dep,eps);

% land mask
mask = zeros(ni,nj);
aa = find(ht > 0);
mask(aa) = 1;

%define a patch to advect
TF = zeros(ni,nj);
TF2 = zeros(ni,nj);

%Natasha NAtl region
TF(206:295,150:177) = 100.0;

TF = TF .* mask;
TF0 = TF;
total_mass(1) = sum(TF(:).*area(:));

%% Following Advect_upwind_2D
dt = 60*60*(1);
ntime = 365 * (60*60*24) / dt;

%I pulled out one more on each side of Natasha's region
iids = [205:296];
jids = [149:178];
uvel = zeros(ni,nj);
vvel = zeros(ni,nj);
uvel(iids,jids) = uh50(:,:,41);
vvel(iids,jids) = vh50(:,:,41);
K = 600.0;

fe = zeros(ni,nj);
fn = zeros(ni,nj);
dfe = zeros(ni,nj);
dfn = zeros(ni,nj);
gradTi = zeros(ni,nj);
gradTj = zeros(ni,nj);
upwind = zeros(ni,nj);
dupwind = zeros(ni,nj);

%% Advection loop
for n = 1:ntime
    n
    
    % Calculate biomass gradient
    %Gradient i
    for j=jsd:jed
        for i=isd:ied
            if (i == ied)
                gradTi(i,j) = (TF(isd,j) - TF(i,j)) ./ dyte(i,j) *mask(i,j)*mask(isd,j);
            else
                gradTi(i,j) = (TF(i+1,j) - TF(i,j)) ./ dyte(i,j) *mask(i,j)*mask(i+1,j);
            end
        end
    end
    %Gradient j
    for j=jsd:jed
        for i=isd:ied
            if (j < jed)
                gradTj(i,j) = (TF(i,j+1) - TF(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(i,j+1);
            else
                gradTj(i,j) = (TF(ni-i+1,j) - TF(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(ni-i+1,j);
            end
        end
    end
    gradT = (gradTi + gradTj); %.* mask;
    
    diffusiv = 0.5*K;
    kpos     = diffusiv + abs(diffusiv);
    kneg     = diffusiv - abs(diffusiv);
    
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
                    dfe(i,j)  = dyte(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(isd,j)) ...
                        .* mask(i,j) .* mask(isd,j);
                else
                    fe(i,j) = dyte(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(i+1,j)/dep(i+1,j))* ...
                        mask(i,j)*mask(i+1,j);
                    dfe(i,j)  = dyte(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(i+1,j)) ...
                        .* mask(i,j) .* mask(i+1,j);
                end
                
            end
        end
    end
    
    % Northward flux
    for j = jsd:jed
        for i = isd:ied
            velocity = 0.5*vvel(i,j);
            upos = velocity + abs(velocity);
            uneg = velocity - abs(velocity);
            
            % define only for ocean cells
            if (mask(i,j) > 0)
                
                if (j < jed)
                    fn(i,j) = dxtn(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(i,j+1)/dep(i,j+1))* ...
                        mask(i,j)*mask(i,j+1);
                    dfn(i,j)  = dxtn(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(i,j+1)) .* ...
                        mask(i,j) .* mask(i,j+1);
                else
                    fn(i,j) = dxtn(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(ni-i+1,j)/dep(ni-i+1,j))* ...
                        mask(i,j)*mask(ni-i+1,j);
                    dfn(i,j)  = dxtn(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(ni-i+1,j)) .* ...
                        mask(i,j) .* mask(ni-i+1,j);
                end
                
            end
        end
    end
    
    % Combine fluxes
    for j = jsd:jed
        for i = isd:ied
            if (j > 1)
                if (i > 1)
                    upwind(i,j) = mask(i,j).*(fe(i-1,j)-fe(i,j)+fn(i,j-1)-fn(i,j));
                    dupwind(i,j) = mask(i,j).*(dfe(i-1,j)-dfe(i,j)+dfn(i,j-1)-dfn(i,j));
                else
                    upwind(i,j) = mask(i,j).*(fe(ied,j)-fe(i,j)+fn(i,j-1)-fn(i,j));
                    dupwind(i,j) = mask(i,j).*(dfe(ied,j)-dfe(i,j)+dfn(i,j-1)-dfn(i,j));
                end
            end
        end
    end
    
    % Update tracers
    for j = jsd:jed
        for i = isd:ied
            TF2(i,j) = TF(i,j) + (dt*upwind(i,j))/area(i,j) - (dt*dupwind(i,j))/area(i,j);
        end
    end
    total_mass(n+1) = sum(TF2(:).*area(:));
    
    % Plot, do mass balance, reset tracer fields
    if n == 1
        figure(1)
        surf(geolon_t,geolat_t,TF); view(2); shading interp; caxis([0 100]);
        print('-dpng',[fpath 'POEM_adv_diff_test_tracer_' cname '_day0.png'])
        
        figure(2)
        m_proj('stereographic','lat',90,'long',30,'radius',30);
        m_pcolor(geolon_t,geolat_t,TF);
        shading interp
        colorbar
        %colormap('jet')
        caxis([0 100])
        m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        title('Tracer day 1')
        print('-dpng',[fpath 'POEM_adv_diff_test_tracer_arcticproj_' cname '_day0.png'])
    end
    
    if n == ntime
        figure(13)
        clf
        surf(geolon_t,geolat_t,TF2); view(2); shading interp;
        caxis([0 100]); %caxis([0 2e3]) to see how big instabilities are
        pdiff = 100*(total_mass(end) - total_mass(1))/total_mass(1);
        title(['%diff = ', num2str(pdiff,'%10.3e')]);
        print('-dpng',[fpath 'POEM_adv_diff_test_tracer_' cname '_day365.png'])
    end
    
    TF = TF2;
    
end

% Arctic projection of tracer
figure(14)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,TF2);
shading interp
colorbar
%colormap('jet')
caxis([0 100])
m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
title('Tracer day 365')
print('-dpng',[fpath 'POEM_adv_diff_test_tracer_arcticproj_' cname '_day365.png'])

save(['/Volumes/GFDL/CSV/advect_tests/POEM_adv_diff_test_' cname '.mat'],'TF0','TF2','pdiff','total_mass')

