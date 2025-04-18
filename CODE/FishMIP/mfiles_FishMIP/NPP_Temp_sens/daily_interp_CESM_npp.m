% Make mat files of interpolated time series from CESM

%NPP control: All forcings pre-industrial (1860-2100) except for temperature. 
%For temperature use historical and rcp85 forcings one after the other 
%(run historical 1860-2005, then RCP8.5 2006-2100)

clear all
close all

fpath='/Volumes/GFDL/Fish-MIP/CESM/PreIndust/';
hpath='/Volumes/GFDL/Fish-MIP/CESM/Hist/';
rpath='/Volumes/GFDL/Fish-MIP/CESM/RCP85/';
spath='/Volumes/GFDL/Fish-MIP/CESM/NPP_cont/';

tstart = 185001:1000:210001;
tend = 185912:1000:210912;
tend(end)=210012;

%each file has 10 years of 12 months, except last is 1 year of 12 months
mos = 10*12*(length(tstart) - 1) + 12;
mstart = 1:12:mos;
mend = 12:12:mos;

yrs= 1850:2100;
Tdays=1:365;
Time=Tdays(15:30:end);

%% Units
%poc flux: mmol C/m^3 cm/s
%szoo: mmol C/m^3
%lzoo: mmol C/m^3
%tp: degC
%tb: degC

%% Altered variables - temp
load([hpath 'cesm_hist_tp_100_monthly_185001-200512.mat'],'tp_100');
load([hpath 'cesm_hist_tbtm_monthly_185001-200512.mat'],'tbtm');
hist_tp_100 = tp_100;
hist_tbtm   = tbtm;
clear tp_100 tbtm

load([rpath 'cesm_rcp85_tp_100_monthly_200601-210012.mat'],'tp_100');
load([rpath 'cesm_rcp85_tbtm_monthly_200601-210012.mat'],'tbtm');
rcp_tp_100 = tp_100;
rcp_tbtm   = tbtm;
clear tp_100 tbtm

tp_100 = cat(3, hist_tp_100, rcp_tp_100);
tbtm   = cat(3, hist_tbtm, rcp_tbtm);

clear hist_tp_100 hist_tbtm rcp_tp_100 rcp_tbtm

%% Control variables - all related to NPP
load([fpath 'cesm_pi_szoo_100_monthly_',num2str(tstart(1)),'-',...
    num2str(tend(end)),'.mat'],'nsmz_100');
load([fpath 'cesm_pi_lzoo_100_monthly_',num2str(tstart(1)),'-',...
    num2str(tend(end)),'.mat'],'nlgz_100');
load([fpath 'cesm_pi_det_btm_monthly_',num2str(tstart(1)),'-',...
        num2str(tend(end)),'.mat'],'poc_btm');

%%
for y = 1:length(yrs)
    yr = yrs(y)
    
    Tp = tp_100(:,:,mstart(y):mend(y));
    Tb = tbtm(:,:,mstart(y):mend(y));
    Zm = nsmz_100(:,:,mstart(y):mend(y));
    Zl = nlgz_100(:,:,mstart(y):mend(y));
    det= poc_btm(:,:,mstart(y):mend(y));
    
    % index of water cells
    [ni,nj,nt] = size(Tp);
    WID = find(~isnan(Tp(:,:,1)));  % spatial index of water cells
    NID = length(WID);              % number of water cells
    
    % setup FEISTY data files
    D_Tp  = zeros(NID,365);
    D_Tb  = zeros(NID,365);
    D_Zm  = zeros(NID,365);
    D_Zl  = zeros(NID,365);
    D_det = zeros(NID,365);
    
    %% interpolate to daily resolution
    for j = 1:NID
        % indexes
        [m,n] = ind2sub([ni,nj],WID(j)); % spatial index of water cell
        
        % pelagic temperature (in Celcius)
        Y = squeeze(Tp(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_Tp(j,:) = yi;
        
        % bottom temperature (in Celcius)
        Y = squeeze(Tb(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_Tb(j,:) = yi;
        
        % medium zoo: from mmolC m-3 to g(WW) m-2
        % 1e-3 mol in 1 mmol
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        % mult by 10 m depth interval for m-3 to m-2
        Y = squeeze(Zm(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_Zm(j,:) = yi * 1e-3 * 12.01 * 9.0 * 10;
        
        % large zoo: from mmolC m-3 to g(WW) m-2
        Y = squeeze(Zl(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_Zl(j,:) = yi * 1e-3  * 12.01 * 9.0 * 10;
        
        % detrital flux to benthos: from mmolC m-3 cm s-1 to g(WW) m-2 d-1
        %poc flux: mmol C/m^3 cm/s
        % 1e-3 mol in 1 mmol
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        % 1e-2 m in cm
        % 60*60*24 sec in a day
        Y = squeeze(det(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_det(j,:) = yi * 1e-3  * 12.01 * 9.0 * 1e-2 * 60 * 60 * 24;
        
    end
    
    % Negative biomass or mortality loss from interp
    D_Zm(D_Zm<0) = 0.0;
    D_Zl(D_Zl<0) = 0.0;
    D_det(D_det<0) = 0.0;
    
    CESM.Tp = D_Tp;
    CESM.Tb = D_Tb;
    CESM.Zm = D_Zm;
    CESM.Zl = D_Zl;
    CESM.det = D_det;
    
    % save
    save([spath 'Data_cesm_npp_cont_' num2str(yr) '.mat'], 'CESM');
    
end
