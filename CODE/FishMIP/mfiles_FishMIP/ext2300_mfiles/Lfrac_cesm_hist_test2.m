% Make mat files of interpolated time series from CESM
% Hist 1950-2014
% 200m vertical integrations

clear 
close all

fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/hist/';
ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/2300/testing_forcing/';

%% Units
%poc flux: mol C m-2 s-1
%zoo: mol C m-2
%tp: degC
%tb: degC

load([fpath 'cesm2_hist_det_monthly_1850_2014.mat'],'det');
load([fpath 'cesm2_hist_phyc_150_monthly_1850_2014.mat'],'phyc_150');
load([fpath 'cesm2_hist_diat_150_monthly_1850_2014.mat'],'diat_150');
load([fpath 'cesm2_hist_zooc_150_monthly_1850_2014.mat']);


%%
zooc_150(zooc_150 > 1.0e19) = nan;
phyc_150(phyc_150 > 1.0e19) = nan;
diat_150(diat_150 > 1.0e19) = nan;
det(det > 1.0e19) = nan;

%% Calc zmeso from diat frac
Lfrac = double(diat_150) ./ double(phyc_150);
Lfrac(Lfrac>1) = 1.0;
Lfrac(Lfrac<0) = 0.0;

zmeso_150_v1 = Lfrac .* double(zooc_150);

zmeso_150_v2 = (Lfrac .* double(zooc_150)) + (0.30 * (1-Lfrac) .* double(zooc_150)); %IPSL pref

zmeso_150_v3 = (Lfrac .* double(zooc_150)) + (0.4286 * (1-Lfrac) .* double(zooc_150)); %UKESM pref

%%
[ni,nj,nt] = size(det);
mos = nt;
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

yrs = floor(yr(1)):floor(yr(end));
Tdays=1:365;

% yrs = 1850:2014;
% Time=Tdays(15:30:end);

%% index of water cells
test4 = squeeze(double(det(:,:,70)));

% Use btm det to find ocean cells
WID = find(~isnan(test4(:)));
NID = length(WID);


[ni,nj] = size(test4);

%%
for y = 1
    ytime = yrs(y);

    if y==1
        range = mstart(y):(mend(y)+1);
        Time=15:30:395;
    elseif y==nyrs
        range = (mstart(y)-1):mend(y);
        Time=-15:30:365;
    else
        range = (mstart(y)-1):(mend(y)+1);
        Time=-15:30:395;
    end

    Zm1 = double(zmeso_150_v1(:,:,range));
    Zm2 = double(zmeso_150_v2(:,:,range));
    Zm3 = double(zmeso_150_v3(:,:,range));

    % setup FEISTY data files
    D_Zm1 = nan*zeros(NID,365);
    D_Zm2 = nan*zeros(NID,365);
    D_Zm3 = nan*zeros(NID,365);

    %% interpolate to daily resolution
    for j = 1:NID
        % indexes
        [m,n] = ind2sub([ni,nj],WID(j)); % spatial index of water cell

        % meso zoo: from molC m-2 to g(WW) m-2
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(Zm1(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Zm1(j,:) = yi * 12.01 * 9.0;

        
        Y = squeeze(Zm2(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Zm2(j,:) = yi * 12.01 * 9.0;

        Y = squeeze(Zm3(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Zm3(j,:) = yi * 12.01 * 9.0;

    end

    % Negative biomass or mortality loss from interp
    D_Zm1(D_Zm1<0) = 0.0;
    D_Zm2(D_Zm2<0) = 0.0;
    D_Zm3(D_Zm3<0) = 0.0;

    ESM.Zm1 = D_Zm1;
    ESM.Zm2 = D_Zm2;
    ESM.Zm3 = D_Zm3;

end


%% Means over all grid cells
[ni,nj,nt] = size(det);
Zm1 = double(reshape(zmeso_150_v1,ni*nj,nt));
Zm2 = double(reshape(zmeso_150_v2,ni*nj,nt));
Zm3 = double(reshape(zmeso_150_v3,ni*nj,nt));

Zm1 = Zm1(WID,:);
Zm2 = Zm2(WID,:);
Zm3 = Zm3(WID,:);

hist_Zm1 = mean(Zm1);
hist_Zm2 = mean(Zm2);
hist_Zm3 = mean(Zm3);

%%
figure
plot(yr,hist_Zm1,'color',[0.75 0 0.5]); hold on
plot(yr,hist_Zm2,'color',[0 0.5 0.75]); hold on
plot(yr,hist_Zm3,'k')

%% save means
hist_yr = yr;
save([fpath 'Lfrac_zoo_test2_cesm2_hist_monthly_1850.mat'],...
    'hist_Zm1','hist_Zm2','hist_Zm3','hist_yr','ESM');

zESM = ESM;

%% Map data

cpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/';
hpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/hist/';

load([cpath 'gridspec_cesm2_cmip6_2300.mat']);
load([cpath 'Data_grid_cesm2_cmip6_2300.mat']);
load([hpath 'Data_cesm_hist_daily_1850.mat'])

[ni,nj]=size(LON);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
clatlim=[plotminlat plotmaxlat];
clonlim=[plotminlon plotmaxlon];

load coastlines

%%
zESM.dZm1 = 10 .^ (-2.925 + 1.964.*log10(zESM.Zm1+eps) + 1.958e-2.*ESM.Tp);
zESM.dZm2 = real(10 .^ (-2.925 + 1.964.*log10(zESM.Zm2+eps) + 1.958e-2.*ESM.Tp));
zESM.dZm3 = real(10 .^ (-2.925 + 1.964.*log10(zESM.Zm3+eps) + 1.958e-2.*ESM.Tp));

c_Z1 = mean(zESM.Zm1,2);
c_Z2 = mean(zESM.Zm2,2);
c_Z3 = mean(zESM.Zm3,2);

c_Zl1 = mean(zESM.dZm1,2);
c_Zl2 = mean(zESM.dZm2,2);
c_Zl3 = mean(zESM.dZm3,2);

[mi,mj]=size(LON);
cZ1=NaN*ones(mi,mj);
cZ2=NaN*ones(mi,mj);
cZ3=NaN*ones(mi,mj);

cZl1=NaN*ones(mi,mj);
cZl2=NaN*ones(mi,mj);
cZl3=NaN*ones(mi,mj);


cZ1(GRD.ID) =c_Z1;
cZ2(GRD.ID) =c_Z2;
cZ3(GRD.ID) =c_Z3;

cZl1(GRD.ID) =c_Zl1;
cZl2(GRD.ID) =c_Zl2;
cZl3(GRD.ID) =c_Zl3;

%%
mod = 'CESM';
f1 = figure('Units','inches','Position',[1 3 5 5]);

subplot('Position',[0.05 0.68 0.45 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(cZ1))
cmocean('tempo')
clim([0 2])
colorbar
title([mod ' Zmeso1 (gWW/m2)'])
%patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.68 0.45 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(cZl1))
cmocean('tempo')
clim([-2 1])
colorbar
title([mod ' Zmeso1 loss (gWW/m2/d)'])
%h2=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.05 0.37 0.45 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(cZ2))
cmocean('tempo')
clim([0 2])
colorbar
title([mod ' Zmeso2 (gWW/m2)'])
%h3=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.37 0.45 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(cZl2))
cmocean('tempo')
clim([-2 1])
colorbar
title([mod ' Zmeso2 loss (gWW/m2/d)'])
%h4=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.05 0.06 0.45 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(cZ3))
cmocean('tempo')
clim([0 2])
colorbar
title([mod ' Zmeso3 (gWW/m2/d)'])
%h5=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.06 0.45 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(cZl3))
cmocean('tempo')
clim([-2 1])
colorbar
title([mod ' Zmeso3 loss (gWW/m2/d)'])
%h5=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'Map_CESM2-WACCM_1850_Lfrac_zoop_test2.png'])

% More mesozoo biomass, but loss rates don't look very different

%%
figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(cZl2-cZl1))
cmocean('tempo')
clim([0 0.02])
colorbar
title([mod ' Zmeso Loss diff (gWW/m2)'])


%% Plot t.s. with other ESMs

Chist_yr = movmean(hist_yr,12);
Chist_Zm1 = movmean(hist_Zm1,12);
Chist_Zm2 = movmean(hist_Zm2,12);
Chist_Zm3 = movmean(hist_Zm3,12);

clear hist_yr

%% Hist UKESM
hpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';

load([hpath 'Means_ukesm_hist_monthly.mat'],'hist_yr','hist_Zm');

Uhist_yr = movmean(hist_yr,12);
Uhist_Zm = movmean(hist_Zm,12);

clear hist_yr hist_Zm

%% Hist IPSL
hpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/hist/';

load([hpath 'Means_ipsl_hist_monthly.mat'],'hist_yr','hist_Zm');

Ihist_yr = movmean(hist_yr,12);
Ihist_Zm = movmean(hist_Zm,12);

clear hist_yr hist_Zm

%%
figure
plot(Uhist_yr,Uhist_Zm,'b','LineWidth',2); hold on
plot(Ihist_yr,Ihist_Zm,'r','LineWidth',2); hold on
plot(Chist_yr,Chist_Zm1,'color',[0 0.75 0.5],'LineWidth',2);
plot(Chist_yr,Chist_Zm2,'color',[0.75 0 0.5],'LineWidth',2);
plot(Chist_yr,Chist_Zm3,'color',[0 0.25 0.75],'LineWidth',2);
title('Mesozoo')
legend('UKESM','IPSL','CESM0','CESM30','CESM43')
legend('location','northwest')
print('-dpng',[ppath 'ts_ESMs_hist_Lfrac_zoop_test2.png'])














