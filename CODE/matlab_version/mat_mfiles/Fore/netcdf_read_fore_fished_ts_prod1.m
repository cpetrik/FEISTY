% POEM output at all locations

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];

%nt=12*95;

%% SP
ncid = netcdf.open([fpath 'Forecast_',harv,'_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(biomass);
%%
SP.bio = biomass;
SP.prod = prod;
% SP.rec = rec;
% SP.clev = clev;
% SP.con = con;
% SP.die = die;
% SP.gamma = gamma;
% SP.nu = nu;

clear biomass clev con die gamma nu rec time prod

% SF
ncid = netcdf.open([fpath 'Forecast_',harv,'_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
SF.prod = prod(:,1:nt);
% SF.rec = rec(:,1:nt);
% SF.clev = clev(:,1:nt);
% SF.con = con(:,1:nt);
% SF.die = die(:,1:nt);
% SF.gamma = gamma(:,1:nt);
% SF.nu = nu(:,1:nt);

clear biomass clev con die gamma nu rec time prod

% SD
ncid = netcdf.open([fpath 'Forecast_',harv,'_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
SD.prod = prod;
% SD.rec = rec;
% SD.clev = clev;
% SD.con = con;
% SD.die = die;
% SD.gamma = gamma;
% SD.nu = nu;

clear biomass clev con die gamma nu rec time prod

% MP
ncid = netcdf.open([fpath 'Forecast_',harv,'_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
MP.prod = prod;
MP.yield = yield;
% MP.rec = rec;
% MP.clev = clev;
% MP.con = con;
% MP.die = die;
% MP.gamma = gamma;
% MP.nu = nu;

clear biomass clev con die gamma nu rec time prod yield

% MF
ncid = netcdf.open([fpath 'Forecast_',harv,'_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
MF.prod = prod;
MF.yield = yield;
% MF.rec = rec;
% MF.rep = rep;
% MF.clev = clev;
% MF.con = con;
% MF.die = die;
% MF.gamma = gamma;
% MF.nu = nu;

clear biomass clev con die gamma nu rec rep time prod yield

% MD
ncid = netcdf.open([fpath 'Forecast_',harv,'_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
MD.prod = prod;
MD.yield = yield;
% MD.rec = rec;
% MD.clev = clev;
% MD.con = con;
% MD.die = die;
% MD.gamma = gamma;
% MD.nu = nu;

clear biomass clev con die gamma nu rec time prod yield

% LP
ncid = netcdf.open([fpath 'Forecast_',harv,'_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
LP.prod = prod;
LP.yield = yield;
% LP.rec = rec;
% LP.rep = rep;
% LP.clev = clev;
% LP.con = con;
% LP.die = die;
% LP.gamma = gamma;
% LP.nu = nu;

clear biomass clev con die gamma nu rec rep time prod yield

% LD
ncid = netcdf.open([fpath 'Forecast_',harv,'_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
LD.prod = prod;
LD.yield = yield;
% LD.rec = rec;
% LD.rep = rep;
% LD.clev = clev;
% LD.con = con;
% LD.die = die;
% LD.gamma = gamma;
% LD.nu = nu;

clear biomass clev con die gamma nu rec rep time prod yield

% Benthic material
ncid = netcdf.open([fpath 'Forecast_',harv,'_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
clear biomass

%% Take means
% Every year
st=1:12:length(time);
en=12:12:length(time);

for n=1:length(st)
    sp_prod1(:,n)=nanmean(SP.prod(:,st(n):en(n)),2);
    sf_prod1(:,n)=nanmean(SF.prod(:,st(n):en(n)),2);
    sd_prod1(:,n)=nanmean(SD.prod(:,st(n):en(n)),2);
    mp_prod1(:,n)=nanmean(MP.prod(:,st(n):en(n)),2);
    mf_prod1(:,n)=nanmean(MF.prod(:,st(n):en(n)),2);
    md_prod1(:,n)=nanmean(MD.prod(:,st(n):en(n)),2);
    lp_prod1(:,n)=nanmean(LP.prod(:,st(n):en(n)),2);
    ld_prod1(:,n)=nanmean(LD.prod(:,st(n):en(n)),2);
end

%% Save means
save([fpath 'Means_fore_',harv,'_' cfile '.mat'],...
'sf_prod1','sp_prod1','sd_prod1',...
'mf_prod1','mp_prod1','md_prod1',...
'lp_prod1','ld_prod1','-append');
