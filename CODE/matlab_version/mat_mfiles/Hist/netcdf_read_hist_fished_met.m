% FEISTY output at all locations
% ESM2M Historic simulation
% Output metabolism

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/Historic_ESM2M/'];

%% SP
ncid = netcdf.open([fpath 'Historic_' harv '_met_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(met);

SP.met = met;
clear met

%% SF
ncid = netcdf.open([fpath 'Historic_' harv '_met_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

SF.met = met(:,1:nt);
clear met

%% SD
ncid = netcdf.open([fpath 'Historic_' harv '_met_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

SD.met = met;
clear met

% MP
ncid = netcdf.open([fpath 'Historic_' harv '_met_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

MP.met = met;
clear met

% MF
ncid = netcdf.open([fpath 'Historic_' harv '_met_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

MF.met = met;
clear met

% MD
ncid = netcdf.open([fpath 'Historic_' harv '_met_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

MD.met = met;
clear met

% LP
ncid = netcdf.open([fpath 'Historic_' harv '_met_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

LP.met = met;
clear met

% LD
ncid = netcdf.open([fpath 'Historic_' harv '_met_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

LD.met = met;
clear met


%% Take means and totals
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%Time
sp_tmet=mean(SP.met,1);
sf_tmet=mean(SF.met,1);
sd_tmet=mean(SD.met,1);
mp_tmet=mean(MP.met,1);
mf_tmet=mean(MF.met,1);
md_tmet=mean(MD.met,1);
lp_tmet=mean(LP.met,1);
ld_tmet=mean(LD.met,1);

%% 50 yrs (1951-2000)
y = 1860+(1/12):(1/12):2005;
yr50=find(y>=1951 & y<2001);
sp_met50=nanmean(SP.met(:,yr50),2);
sf_met50=nanmean(SF.met(:,yr50),2);
sd_met50=nanmean(SD.met(:,yr50),2);
mp_met50=nanmean(MP.met(:,yr50),2);
mf_met50=nanmean(MF.met(:,yr50),2);
md_met50=nanmean(MD.met(:,yr50),2);
lp_met50=nanmean(LP.met(:,yr50),2);
ld_met50=nanmean(LD.met(:,yr50),2);

%% Seasonal climatology
nmo = nt;
nyr = nt/12;

all_clim = nan*ones(ni,12);
sf_clim = all_clim;
sp_clim = all_clim;
sd_clim = all_clim;
mf_clim = all_clim;
mp_clim = all_clim;
md_clim = all_clim;
lp_clim = all_clim;
ld_clim = all_clim;
for m = 1:12
    mo = m:12:nmo;
    sf_clim(:,m) = nanmean(SF.met(:,mo),2);
    sp_clim(:,m) = nanmean(SP.met(:,mo),2);
    sd_clim(:,m) = nanmean(SD.met(:,mo),2);
    mf_clim(:,m) = nanmean(MF.met(:,mo),2);
    mp_clim(:,m) = nanmean(MP.met(:,mo),2);
    md_clim(:,m) = nanmean(MD.met(:,mo),2);
    lp_clim(:,m) = nanmean(LP.met(:,mo),2);
    ld_clim(:,m) = nanmean(LD.met(:,mo),2);
end

%% quick plot
figure
subplot(2,2,1)
plot(1:12,mf_clim(1,:))
subplot(2,2,2)
plot(1:12,mf_clim(100,:))
subplot(2,2,3)
plot(1:12,mf_clim(1000,:))
subplot(2,2,4)
plot(1:12,mf_clim(10000,:))

figure
subplot(2,2,1)
plot(1:12,mf_clim(100,:))
subplot(2,2,2)
plot(1:12,lp_clim(100,:))
subplot(2,2,3)
plot(1:12,ld_clim(100,:))
subplot(2,2,4)
plot(1:12,sf_clim(100,:))

%%
save([fpath 'Means_Historic_' harv '_met_' cfile '.mat'],...
    'y','yr50',...
    'sf_tmet','sp_tmet','sd_tmet',...
    'mf_tmet','mp_tmet','md_tmet',...
    'lp_tmet','ld_tmet',...
    'sf_met50','sp_met50','sd_met50',...
    'mf_met50','mp_met50','md_met50',...
    'lp_met50','ld_met50',...
    'sf_clim','sp_clim','sd_clim',...
    'mf_clim','mp_clim','md_clim',...
    'lp_clim','ld_clim');
