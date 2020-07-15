% FEISTY output at all locations
% Fraction of zoop hp loss consumed

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];

% MZ
ncid = netcdf.open([fpath 'Historic_' harv '_mzoo.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);
%%
MZ.frac = fraction;
clear fraction time

% LZ
ncid = netcdf.open([fpath 'Historic_' harv '_lzoo.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Historic_' harv '_bfrac.nc'],'NC_NOWRITE');
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

%total times it happens
md_ttfrac=nansum(MD.catch,1);
lp_ttfrac=nansum(LP.catch,1);
ld_ttfrac=nansum(LD.catch,1);

%% 50 yrs (1951-2000)
y = 1860+(1/12):(1/12):2005;
yr50=find(y>=1951 & y<2001);
mf_mc50=nanmean(MF.catch(:,yr50),2);
mp_mc50=nanmean(MP.catch(:,yr50),2);
md_mc50=nanmean(MD.catch(:,yr50),2);
lp_mc50=nanmean(LP.catch(:,yr50),2);
ld_mc50=nanmean(LD.catch(:,yr50),2);
mf_tc50=nansum(MF.catch(:,yr50),2);
mp_tc50=nansum(MP.catch(:,yr50),2);
md_tc50=nansum(MD.catch(:,yr50),2);
lp_tc50=nansum(LP.catch(:,yr50),2);
ld_tc50=nansum(LD.catch(:,yr50),2);

% 1990-1995 (Climatology)
lyr=find(y>=1990 & y<1995);
%Means
mf_mc5=nanmean(MF.catch(:,lyr),2);
mp_mc5=nanmean(MP.catch(:,lyr),2);
md_mc5=nanmean(MD.catch(:,lyr),2);
lp_mc5=nanmean(LP.catch(:,lyr),2);
ld_mc5=nanmean(LD.catch(:,lyr),2);
mf_tc5=nansum(MF.catch(:,lyr),2);
mp_tc5=nansum(MP.catch(:,lyr),2);
md_tc5=nansum(MD.catch(:,lyr),2);
lp_tc5=nansum(LP.catch(:,lyr),2);
ld_tc5=nansum(LD.catch(:,lyr),2);

% 1990 (Climatology)
yr1=find(y>=1990 & y<1991);
%Totals
mp_tc90=nansum(MP.catch(:,yr1),2);
mf_tc90=nansum(MF.catch(:,yr1),2);
md_tc90=nansum(MD.catch(:,yr1),2);
lp_tc90=nansum(LP.catch(:,yr1),2);
ld_tc90=nansum(LD.catch(:,yr1),2);
mp_mc90=nanmean(MP.catch(:,yr1),2);
mf_mc90=nanmean(MF.catch(:,yr1),2);
md_mc90=nanmean(MD.catch(:,yr1),2);
lp_mc90=nanmean(LP.catch(:,yr1),2);
ld_mc90=nanmean(LD.catch(:,yr1),2);

%% Every 5 years
st=1:60:length(time);
en=60:60:length(time);

for n=1:length(st)

    mp_mc(:,n)=nanmean(MP.catch(:,st(n):en(n)),2);
    mf_mc(:,n)=nanmean(MF.catch(:,st(n):en(n)),2);
    md_mc(:,n)=nanmean(MD.catch(:,st(n):en(n)),2);
    lp_mc(:,n)=nanmean(LP.catch(:,st(n):en(n)),2);
    ld_mc(:,n)=nanmean(LD.catch(:,st(n):en(n)),2);

    mp_tc(:,n)=nansum(MP.catch(:,st(n):en(n)),2);
    mf_tc(:,n)=nansum(MF.catch(:,st(n):en(n)),2);
    md_tc(:,n)=nansum(MD.catch(:,st(n):en(n)),2);
    lp_tc(:,n)=nansum(LP.catch(:,st(n):en(n)),2);
    ld_tc(:,n)=nansum(LD.catch(:,st(n):en(n)),2);
end

% Every year
st=1:12:length(time);
en=12:12:length(time);

mf_tyc = nan*ones(ni,nyr);
mp_tyc = nan*ones(ni,nyr);
md_tyc = nan*ones(ni,nyr);
lp_tyc = nan*ones(ni,nyr);
ld_tyc = nan*ones(ni,nyr);
for n=1:length(st)
    
    mp_tyc(:,n)=nansum(MP.catch(:,st(n):en(n)),2);
    mf_tyc(:,n)=nansum(MF.catch(:,st(n):en(n)),2);
    md_tyc(:,n)=nansum(MD.catch(:,st(n):en(n)),2);
    lp_tyc(:,n)=nansum(LP.catch(:,st(n):en(n)),2);
    ld_tyc(:,n)=nansum(LD.catch(:,st(n):en(n)),2);
end

tmn = mf_tyc + mp_tyc + md_tyc + lp_tyc + ld_tyc;
stmn = sum(tmn);

mp_tsyc = nansum(mp_tyc);
mf_tsyc = nansum(mf_tyc);
md_tsyc = nansum(md_tyc);
lp_tsyc = nansum(lp_tyc);
ld_tsyc = nansum(ld_tyc);

%%
save([fpath 'Means_Historic_' harv '_' cfile '.mat'],'time',...
    'mf_tyc','mp_tyc','md_tyc','lp_tyc','ld_tyc',...
    'mf_tsyc','mp_tsyc','md_tsyc','lp_tsyc','ld_tsyc',...
    'units_yield','units_catch','-append');







