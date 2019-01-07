% FEISTY output CC
% 3 km model

clear all
close all


cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
%harv = 'pristine';
harv = 'All_fish03';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/CalCurr/'];

%% Benthic material
ncid = netcdf.open([fpath 'Historic3km_' harv '_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

nt = length(time);

Bent.bio = biomass;
BENT.bio = biomass(:,nt);
clear biomass

% SP
ncid = netcdf.open([fpath 'Historic3km_' harv '_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.bio = biomass;
Sml_p.bio = biomass(:,nt);
clear biomass 

% SF
ncid = netcdf.open([fpath 'Historic3km_' harv '_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
Sml_f.bio = biomass(:,nt);
clear biomass 

% SD
ncid = netcdf.open([fpath 'Historic3km_' harv '_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
Sml_d.bio = biomass(:,nt);
clear biomass 

% MP
ncid = netcdf.open([fpath 'Historic3km_' harv '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
MP.yield = yield;
Med_p.bio = biomass(:,nt);
clear biomass yield

% MF
ncid = netcdf.open([fpath 'Historic3km_' harv '_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
MF.yield = yield;
Med_f.bio = biomass(:,nt);
clear biomass yield

% MD
ncid = netcdf.open([fpath 'Historic3km_' harv '_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
MD.yield = yield;
Med_d.bio = biomass(:,nt);
clear biomass yield

% LP
ncid = netcdf.open([fpath 'Historic3km_' harv '_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
LP.yield = yield;
Lrg_p.bio = biomass(:,nt);
clear biomass yield

% LD
ncid = netcdf.open([fpath 'Historic3km_' harv '_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
LD.yield = yield;
Lrg_d.bio = biomass(:,nt);
clear biomass yield

%% Take means
% Map data
cpath = '/Volumes/GFDL/NEMURO/3km/';
load([cpath 'gridspec_3km.mat']);
load([cpath 'Data_grid_3km_hist.mat']);

coast = find(GRD.DIST <=200);

%% Time
sp_tmean_coast=mean(SP.bio(coast,:),1);
sf_tmean_coast=mean(SF.bio(coast,:),1);
sd_tmean_coast=mean(SD.bio(coast,:),1);
mp_tmean_coast=mean(MP.bio(coast,:),1);
mf_tmean_coast=mean(MF.bio(coast,:),1);
md_tmean_coast=mean(MD.bio(coast,:),1);
lp_tmean_coast=mean(LP.bio(coast,:),1);
ld_tmean_coast=mean(LD.bio(coast,:),1);
b_tmean_coast=mean(Bent.bio(coast,:),1);

mf_tmy_coast=mean(MF.yield(coast,:),1);
mp_tmy_coast=mean(MP.yield(coast,:),1);
md_tmy_coast=mean(MD.yield(coast,:),1);
lp_tmy_coast=mean(LP.yield(coast,:),1);
ld_tmy_coast=mean(LD.yield(coast,:),1);

%%
save([fpath 'Means_Historic3km_' harv '_' cfile '.mat'],...
    'sf_tmean_coast','sp_tmean_coast','sd_tmean_coast',...
    'mf_tmean_coast','mp_tmean_coast','md_tmean_coast',...
    'b_tmean_coast','lp_tmean_coast','ld_tmean_coast','time',...
    'mf_tmy_coast','mp_tmy_coast','md_tmy_coast',...
    'lp_tmy_coast','ld_tmy_coast','coast','-append');

% Save last year for initializing forecast runs
% save([fpath 'Last_mo_Historic3km_' harv '_' cfile '.mat'],'Sml_f','Sml_p','Sml_d',... 
%     'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')

% end





