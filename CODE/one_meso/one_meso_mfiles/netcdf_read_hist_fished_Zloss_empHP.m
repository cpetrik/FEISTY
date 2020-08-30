% Hist
% FEISTY output at all locations
% Fraction of zoop hp loss consumed

close all
clear all

%cfile = 'Dc_Lam700_enc70-b250_m400-b175-k086_c20-b250_D080_A050_nmort1_BE08_CC80_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];

% MZ
ncid = netcdf.open([fpath 'Historic_1meso_empHPloss_' harv '_mzoo.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MZ.fracL = fractionLoss;


%% Take means and totals

% Time
%mean yield per mo
mz_tmfrac=nanmean(MZ.fracL,1);

%% Space
% Last year
lyr=time((end-12+1):end);
mz_mfrac5=nanmean(MZ.fracL(:,lyr),2);

%% Total times it happens
MZ.over = nan*ones(size(MZ.fracL));

MZ.over(MZ.fracL > 1) = ones;

MZ.over(MZ.fracL <= 1) = zeros;

% Time
mz_ttf=nansum(MZ.over,1);

% Space
mz_mtf5=nansum(MZ.over(:,lyr),2);

%% Every year
[ni,nt] = size(MZ.fracL);
nyr = nt/12;
st=1:12:length(time);
en=12:12:length(time);
mz_mtf = nan*ones(ni,nyr);
for n=1:length(st)
    mz_mtf(:,n)=nansum(MZ.over(:,st(n):en(n)),2);
end

%%
save([fpath 'Means_Historic_1meso_empHPloss_' harv '_' cfile '.mat'],'time',...
    'mz_tmfrac','mz_mfrac5','mz_ttf',...
    'mz_mtf5','mz_mtf','-append');
