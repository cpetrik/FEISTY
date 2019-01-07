% Get mean biomasses of temp, det, zoop

clear all
close all

Cdir = '/Volumes/GFDL/NEMURO/10km/';
MASK = ncread([Cdir 'wc12_avgmo_1988.nc'],'MASK');
yrs= 1988:2010;

% index of water cells
[ni,nj] = size(MASK);
WID = find(MASK(:)>0); % spatial index of water cells
NID = length(WID); % number of water cells

% Mortality
%Large zoo: 0.06 m3/(mmolN*day)
%Predatory zoo: 0.08 m3/(mmolN*day)
% Flux of detritus = PON * sinking vel (mmolN/m2/s)

% All units already converted to C, gWW/m2, or gWW/m2/d

%%
nyrs = length(yrs);
D_Tp  = zeros(NID,nyrs);
D_Tb  = zeros(NID,nyrs);
D_Zm  = zeros(NID,nyrs);
D_Zl  = zeros(NID,nyrs);
D_dZm = zeros(NID,nyrs);
D_dZl = zeros(NID,nyrs);
D_det = zeros(NID,nyrs);
for y = 1:nyrs
    yr = yrs(y);
    
    di = [Cdir '/daily10km/Data_10km_' num2str(yr)];
    load([di '.mat'], 'NEMURO');
    
    D_Tp(:,y)  = nanmean(NEMURO.Tp,2);
    D_Tb(:,y)  = nanmean(NEMURO.Tb,2);
    D_Zm(:,y)  = nanmean(NEMURO.Zm,2);
    D_Zl(:,y)  = nanmean(NEMURO.Zl,2);
    D_dZm(:,y) = nanmean(NEMURO.dZm,2);
    D_dZl(:,y) = nanmean(NEMURO.dZl,2);
    D_det(:,y) = nanmean(NEMURO.det,2);
    
end

%%
ptemp_mean_10km = double(nanmean(D_Tp,2));
btemp_mean_10km = double(nanmean(D_Tb,2));
det_mean_10km = double(nanmean(D_det,2));
mz_mean_10km = double(nanmean(D_Zm,2));
lz_mean_10km = double(nanmean(D_Zl,2));
mzloss_mean_10km = double(nanmean(D_dZm,2));
lzloss_mean_10km = double(nanmean(D_dZl,2));

%%
save([Cdir 'nemuro_10km_means.mat'],'WID','D_Tp','D_Tb','D_Zm','D_Zl','D_dZm',...
    'D_dZl','D_det','ptemp_mean_10km','btemp_mean_10km','det_mean_10km',...
    'mz_mean_10km','lz_mean_10km','mzloss_mean_10km','lzloss_mean_10km');






