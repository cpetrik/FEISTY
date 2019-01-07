% Make mat files of interpolated time series from NEMURO 3km model CC
% 1988-2010

clear all
close all

Cdir = '/Volumes/GFDL/NEMURO/3km/';
MASK = ncread([Cdir 'wc15n_avgmo_1988.nc'],'MASK');
yrs= 1988:2010;
Tdays=1:365;
Time=Tdays(15:30:end);

% Mortality
%Large zoo: 0.06 m3/(mmolN*day)
%Predatory zoo: 0.08 m3/(mmolN*day)

% Flux of detritus = PON * sinking vel (mmolN/m2/s)

%%
for y = 1:length(yrs)
    yr = yrs(y)
    
    Tp  = ncread([Cdir 'wc15n_avgmo_',num2str(yr),'.nc'],'T100M');
    TBT = ncread([Cdir 'wc15n_avgmo_',num2str(yr),'.nc'],'TBOT');
    Tb  = squeeze(TBT(:,:,1,:));
    Zm  = ncread([Cdir 'wc15n_avgmo_',num2str(yr),'.nc'],'ZL100M');
    Zl  = ncread([Cdir 'wc15n_avgmo_',num2str(yr),'.nc'],'ZP100M');
    PON = ncread([Cdir 'wc15n_avgmo_',num2str(yr),'.nc'],'PONFLX');
    det = squeeze(PON(:,:,1,:));
    
    Tp(Tp==-1.0000e+34) = nan;
    Tb(Tb==-1.0000e+34) = nan;
    Zm(Zm==-1.0000e+34) = nan;
    Zl(Zl==-1.0000e+34) = nan;
    det(det==-1.0000e+34) = nan;
    
    dZm = 0.06 * Zm;
    dZl = 0.08 * Zl;
        
    % index of water cells
    [ni,nj] = size(MASK);
    WID = find(MASK(:)>0); % spatial index of water cells
    NID = length(WID); % number of water cells
    
    % setup FEISTY data files
    D_Tp  = zeros(NID,365);
    D_Tb  = zeros(NID,365);
    D_Zm  = zeros(NID,365);
    D_Zl  = zeros(NID,365);
    D_dZm = zeros(NID,365);
    D_dZl = zeros(NID,365);
    D_det = zeros(NID,365);
    % D_u = zeros(NID,365);
    % D_v = zeros(NID,365);
    
    % NaN velocities
    % u[find(u.==minimum(u))] = 0.0
    % v[find(v.==minimum(v))] = 0.0
    
    %% interpolate to daily resolution
    for j = 1:NID
        % indexes
        [m,n] = ind2sub([ni,nj],WID(j)); % spatial index of water cell
        
        %     % v currents in m/s
        %     Y = squeeze(V200(m,n,:));
        %     yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        %     D_v(j,:) = yi;
        %
        %     % u currents from m/s to m/d
        %     Y = squeeze(U200(m,n,:));
        %     yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        %     D_u(j,:) = yi;
        
        % pelagic temperature (in Celcius)
        Y = squeeze(Tp(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_Tp(j,:) = yi;
        
        % bottom temperature (in Celcius)
        Y = squeeze(Tb(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_Tb(j,:) = yi;
        
        % medium zoo: from micro mol N m-2 to g(WW) m-2
        % 1e-3 mol in 1 mmol
        % 106/16 mol C in 1 mol N
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(Zm(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_Zm(j,:) = yi * 1e-3 * (106.0/16.0) * 12.01 * 9.0;
        
        % large zoo: from micro mol N m-2 to g(WW) m-2
        Y = squeeze(Zl(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_Zl(j,:) = yi * 1e-3  * (106.0/16.0) * 12.01 * 9.0;
        
        % medium zoo mortality: from micro mol N m-2 d-1 to g(WW) m-2 d-1
        Y = squeeze(dZm(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_dZm(j,:) = yi * 1e-3  * (106.0/16.0) * 12.01 * 9.0;
        
        % large zoo mortality: from micro mol N m-2 d-1 to g(WW) m-2 d-1
        Y = squeeze(dZl(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_dZl(j,:) = yi * 1e-3  * (106.0/16.0) * 12.01 * 9.0;
        
        % detrital flux to benthos: from micro mol C m-2 s-1 to g(WW) m-2 d-1
        Y = squeeze(det(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_det(j,:) = yi * 1e-3  * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
        
    end
    
    % Negative biomass or mortality loss from interp
    D_Zm(D_Zm<0) = 0.0;
    D_Zl(D_Zl<0) = 0.0;
    D_dZm(D_dZm<0) = 0.0;
    D_dZl(D_dZl<0) = 0.0;
    D_det(D_det<0) = 0.0;
    
    NEMURO.Tp = D_Tp;
    NEMURO.Tb = D_Tb;
    NEMURO.Zm = D_Zm;
    NEMURO.Zl = D_Zl;
    NEMURO.dZm = D_dZm;
    NEMURO.dZl = D_dZl;
    NEMURO.det = D_det;
    % NEMURO.U = D_u;
    % NEMURO.V = D_v;
    
    % save
    
    di = [Cdir '/daily3km/Data_3km_' num2str(yr)];
    save([di '.mat'], 'NEMURO');
end



