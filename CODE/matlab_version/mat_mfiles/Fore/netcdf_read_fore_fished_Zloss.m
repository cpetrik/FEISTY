% FEISTY output at all locations
% Fraction of zoop hp loss consumed

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];

% MZ
ncid = netcdf.open([fpath 'Forecast_' harv '_mzoo.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Forecast_' harv '_lzoo.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Forecast_' harv '_bfrac.nc'],'NC_NOWRITE');
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

%% 50 yrs (1951-2000)
y = 2005+(1/12):(1/12):2100;
yr50=find(y>=2051 & y<2100);
mz_mfrac50=nanmean(MZ.frac(:,yr50),2);
lz_mfrac50=nanmean(LZ.frac(:,yr50),2);
b_mfrac50=nanmean(B.frac(:,yr50),2);

% 1990-1995 (Climatology) Means
lyr=find(y>=1990 & y<1995);
mz_mfrac5=nanmean(MZ.frac(:,lyr),2);
lz_mfrac5=nanmean(LZ.frac(:,lyr),2);
b_mfrac5=nanmean(B.frac(:,lyr),2);

% 1990 (Climatology) Means
yr1=find(y>=1990 & y<1991);
mz_mfrac90=nanmean(MZ.frac(:,yr1),2);
lz_mfrac90=nanmean(LZ.frac(:,yr1),2);
b_mfrac90=nanmean(B.frac(:,yr1),2);

% Every year
[ni,nt] = size(LZ.frac);
nyr = nt/12;
st=1:12:length(time);
en=12:12:length(time);
mz_mfrac = nan*ones(ni,nyr);
lz_mfrac = nan*ones(ni,nyr);
b_mfrac = nan*ones(ni,nyr);
for n=1:length(st)
    mz_mfrac(:,n)=nanmean(MZ.frac(:,st(n):en(n)),2);
    lz_mfrac(:,n)=nanmean(LZ.frac(:,st(n):en(n)),2);
    b_mfrac(:,n)=nanmean(B.frac(:,st(n):en(n)),2);
end


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

%50 yrs
mz_mtf50=nansum(MZ.over(:,yr50),2);
lz_mtf50=nansum(LZ.over(:,yr50),2);
b_mtf50=nansum(B.over(:,yr50),2);

% 1990-1995 (Climatology)
mz_mtf5=nansum(MZ.over(:,lyr),2);
lz_mtf5=nansum(LZ.over(:,lyr),2);
b_mtf5=nansum(B.over(:,lyr),2);

% 1990 (Climatology)
mz_mtf90=nansum(MZ.over(:,yr1),2);
lz_mtf90=nansum(LZ.over(:,yr1),2);
b_mtf90=nanmean(B.over(:,yr1),2);

%% Every year
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
save([fpath 'Forecast_RCP85_ESM2M/Means_fore_' harv '_' cfile '.mat'],'time',...
    'mz_tmfrac','mz_mfrac50','mz_mfrac5','mz_mfrac90','mz_mfrac','mz_ttf',...
    'mz_mtf50','mz_mtf5','mz_mtf90','mz_mtf',...
    'lz_tmfrac','lz_mfrac50','lz_mfrac5','lz_mfrac90','lz_mfrac','lz_ttf',...
    'lz_mtf50','lz_mtf5','lz_mtf90','lz_mtf',...
    'b_tmfrac','b_mfrac50','b_mfrac5','b_mfrac90','b_mfrac','b_ttf',...
    'b_mtf50','b_mtf5','b_mtf90','b_mtf','-append');
