% Climatology
% FEISTY output at all locations
% Fraction of zoop hp loss consumed

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];

% MZ
ncid = netcdf.open([fpath 'Climatol_' harv '_mzoo.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MZ.frac = fraction;
clear fraction time

% LZ
ncid = netcdf.open([fpath 'Climatol_' harv '_lzoo.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LZ.frac = fraction;
clear fraction time

% Bent
ncid = netcdf.open([fpath 'Climatol_' harv '_bfrac.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

B.frac = fraction;
clear fraction

%% Take means and totals

% Time
%mean yield per mo
mz_tmfrac=nanmean(MZ.frac,1);
lz_tmfrac=nanmean(LZ.frac,1);
b_tmfrac=nanmean(B.frac,1);

%% Space
% Last year
lyr=time((end-12+1):end);
mz_mfrac5=nanmean(MZ.frac(:,lyr),2);
lz_mfrac5=nanmean(LZ.frac(:,lyr),2);
b_mfrac5=nanmean(B.frac(:,lyr),2);

%% Total times it happens
MZ.over = nan*ones(size(MZ.frac));
LZ.over = nan*ones(size(MZ.frac));
B.over = nan*ones(size(MZ.frac));

MZ.over(MZ.frac > 1) = ones;
LZ.over(LZ.frac > 1) = ones;
B.over(B.frac > 1) = ones;

MZ.over(MZ.frac <= 1) = zeros;
LZ.over(LZ.frac <= 1) = zeros;
B.over(B.frac <= 1) = zeros;

% Time
mz_ttf=nansum(MZ.over,1);
lz_ttf=nansum(LZ.over,1);
b_ttf=nansum(B.over,1);

% Space
mz_mtf5=nansum(MZ.over(:,lyr),2);
lz_mtf5=nansum(LZ.over(:,lyr),2);
b_mtf5=nansum(B.over(:,lyr),2);

%% Every year
[ni,nt] = size(LZ.frac);
nyr = nt/12;
st=1:12:length(time);
en=12:12:length(time);
mz_mtf = nan*ones(ni,nyr);
lz_mtf = nan*ones(ni,nyr);
b_mtf = nan*ones(ni,nyr);
for n=1:length(st)
    mz_mtf(:,n)=nansum(MZ.over(:,st(n):en(n)),2);
    lz_mtf(:,n)=nansum(LZ.over(:,st(n):en(n)),2);
    b_mtf(:,n)=nansum(B.over(:,st(n):en(n)),2);
end

%%
save([fpath 'Climatology/Means_Climatol_' harv '_' cfile '.mat'],'time',...
    'mz_tmfrac','mz_mfrac5','mz_ttf',...
    'mz_mtf5','mz_mtf',...
    'lz_tmfrac','lz_mfrac5','lz_ttf',...
    'lz_mtf5','lz_mtf',...
    'b_tmfrac','b_mfrac5','b_ttf',...
    'b_mtf5','b_mtf','-append');
