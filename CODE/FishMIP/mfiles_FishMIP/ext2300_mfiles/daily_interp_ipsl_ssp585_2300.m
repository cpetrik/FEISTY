% Make mat files of interpolated time series from IPSL
% SSP 585 2101-2300
% 200 m vertical integrations

clear 
close all

fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp585/';

%% Units
%poc flux: mol C m-2 s-1
%zoo: mol C m-2
%tp: degC
%tb: degC

load([fpath 'ipsl_ssp585_temp_btm_monthly_2101_2300.mat'],'temp_btm');
load([fpath 'ipsl_ssp585_temp_200_monthly_2101_2300.mat'],'temp_200');
load([fpath 'ipsl_ssp585_zmeso_200_monthly_2101_2300.mat'],'zmeso_200');
load([fpath 'ipsl_ssp585_det_monthly_2101_2300.mat']); %,'det'

load('/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/gridspec_ipsl_cmip6.mat','deptho')

temp_200(temp_200 > 1.0e19) = nan;
temp_btm(temp_btm > 1.0e19) = nan;
zmeso_200(zmeso_200 > 1.0e19) = nan;
det(det > 1.0e19) = nan;

%%
mos = length(time);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

yrs = floor(yr(1)):floor(yr(end));
Tdays=1:365;

% yrs = 1850:2014;
% Time=Tdays(15:30:end);

%% test that all same orientation
test1 = squeeze(double(temp_200(:,:,70)));
test2 = squeeze(double(temp_btm(:,:,70)));
test3 = squeeze(double(zmeso_200(:,:,70)));
test4 = squeeze(double(det(:,:,70)));

figure
subplot(2,2,1)
pcolor(test1); shading flat
subplot(2,2,2)
pcolor(test2); shading flat
subplot(2,2,3)
pcolor(test3); shading flat
subplot(2,2,4)
pcolor(test4); shading flat

figure
pcolor(deptho); shading flat

%% flip depth
depth = fliplr(deptho);

figure
pcolor(depth); shading flat

%% index of water cells
%make GRD in another file later
% load('/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/Data_grid_ipsl.mat','GRD');
% WID = GRD.ID;
% NID = GRD.N;

% Use btm det to find ocean cells
WID = find(~isnan(test4(:)));
NID = length(WID);

[ni,nj] = size(test4);

%%
for y = 1:nyrs
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
    
    Tp = double(temp_200(:,:,range));
    Tb = double(temp_btm(:,:,range));
    Zm = double(zmeso_200(:,:,range));
    Det= double(det(:,:,range));
    
    % setup FEISTY data files
    D_Tp  = nan*zeros(NID,365);
    D_Tb  = nan*zeros(NID,365);
    D_Zm  = nan*zeros(NID,365);
    D_det = nan*zeros(NID,365);
    
    %% interpolate to daily resolution
    for j = 1:NID
        % indexes
        [m,n] = ind2sub([ni,nj],WID(j)); % spatial index of water cell
        
        % pelagic temperature (in Celcius)
        Y = squeeze(Tp(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Tp(j,:) = yi;
        
        % bottom temperature (in Celcius)
        Y = squeeze(Tb(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Tb(j,:) = yi;
        
        % meso zoo: from molC m-2 to g(WW) m-2
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(Zm(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Zm(j,:) = yi * 12.01 * 9.0;
        
        % detrital flux to benthos: from molC m-2 s-1 to g(WW) m-2 d-1
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        % 60*60*24 sec in a day
        Y = squeeze(Det(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_det(j,:) = yi * 12.01 * 9.0 * 60 * 60 * 24;
        
    end

    % Negative biomass or mortality loss from interp
    D_Zm(D_Zm<0) = 0.0;
    D_det(D_det<0) = 0.0;
    
    ESM.Tp = D_Tp;
    ESM.Tb = D_Tb;
    ESM.Zm = D_Zm;
    ESM.det = D_det;
    
    % save
    save([fpath 'Data_ipsl_ssp585_daily_',num2str(yr),'.mat'], 'ESM');
    
end


%% Means over all grid cells
nt = length(time);

Tp = double(reshape(temp_200,ni*nj,nt));
Tb = double(reshape(temp_btm,ni*nj,nt));
Zm = double(reshape(zmeso_200,ni*nj,nt));
Det= double(reshape(det,ni*nj,nt));

Tp = Tp(WID,:);
Tb = Tb(WID,:);
Zm = Zm(WID,:);
Det= Det(WID,:);

ssp585_Tp2 = mean(Tp);
ssp585_Tb2 = mean(Tb);
ssp585_Zm2 = mean(Zm);
ssp585_Det2 = mean(Det);

%%
figure
subplot(2,2,1)
plot(yr,ssp585_Tp2,'r')

subplot(2,2,2)
plot(yr,ssp585_Tb2,'b')

subplot(2,2,3)
plot(yr,ssp585_Zm2,'color',[0.75 0 0.5])

subplot(2,2,4)
plot(yr,ssp585_Det2,'color',[0 0.5 0.75])

%% save means
ssp585_yr2 = yr;
save([fpath 'Means_ipsl_ssp585_monthly_2101_2300.mat'], 'ssp585_Tp2','ssp585_Tb2',...
    'ssp585_Zm2','ssp585_Det2','ssp585_yr2');

