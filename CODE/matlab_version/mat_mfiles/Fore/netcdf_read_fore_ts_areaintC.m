% POEM output at all locations
% Save area-integrated fractions

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
epath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat','AREA_OCN');
grid = csvread([cpath 'grid_csv.csv']);
area = AREA_OCN(grid(:,1));
area_mat = repmat(area,1,95*12);

%nt=12*95;

%% read netcdfs
% SP
ncid = netcdf.open([fpath 'Forecast_',harv,'_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(biomass);

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

%% Take means and totals

SF.prod(SF.prod(:)<0)=0;
SP.prod(SP.prod(:)<0)=0;
SD.prod(SD.prod(:)<0)=0;
MF.prod(MF.prod(:)<0)=0;
MP.prod(MP.prod(:)<0)=0;
MD.prod(MD.prod(:)<0)=0;
LP.prod(LP.prod(:)<0)=0;
LD.prod(LD.prod(:)<0)=0;

%Time
sp_tamean=nansum(SP.prod.*area,1);
sf_tamean=nansum(SF.prod.*area,1);
sd_tamean=nansum(SD.prod.*area,1);
mp_tamean=nansum(MP.prod.*area,1);
mf_tamean=nansum(MF.prod.*area,1);
md_tamean=nansum(MD.prod.*area,1);
lp_tamean=nansum(LP.prod.*area,1);
ld_tamean=nansum(LD.prod.*area,1);
b_tamean=nansum(Bent.bio.*area,1);

F = SF.prod + MF.prod;
P = SP.prod + MP.prod + LP.prod;
D = SD.prod + MD.prod + LD.prod;
All = F+P+D;

tPD = nansum(P.*area) ./ nansum((P.*area)+(D.*area));
tFD = nansum(F.*area) ./ nansum((F.*area)+(D.*area));
tPelD = nansum((P.*area)+(F.*area)) ./ ...
    nansum((P.*area)+(F.*area)+(D.*area));
tDPel = nansum(D.*area) ./ nansum((P.*area)+(F.*area)+(D.*area));
tDP = nansum(D.*area) ./ nansum((D.*area)+(P.*area));
tDF = nansum(D.*area) ./ nansum((D.*area)+(F.*area));
tFP = nansum(F.*area) ./ nansum((F.*area)+(P.*area));
tPF = nansum(P.*area) ./ nansum((P.*area)+(F.*area));

tP = nansum(P.*area) ./ nansum(All.*area);
tF = nansum(F.*area) ./ nansum(All.*area);
tPel = nansum((P.*area)+(F.*area)) ./ nansum(All.*area);

%%
save([fpath 'Means_fore_',harv,'_' cfile '.mat'],...
    'sf_tamean','sp_tamean','sd_tamean',...
    'mf_tamean','mp_tamean','md_tamean',...
    'lp_tamean','ld_tamean','b_tamean',...
    'tPD','tFD','tPelD','tDPel','tDP','tDF','tFP','tPF',...
    'tP','tF','tPel','-append');

save([epath 'Means_fore_',harv,'_' cfile '.mat'],...
    'sf_tamean','sp_tamean','sd_tamean',...
    'mf_tamean','mp_tamean','md_tamean',...
    'lp_tamean','ld_tamean','b_tamean',...
    'tPD','tFD','tPelD','tDPel','tDP','tDF','tFP','tPF',...
    'tP','tF','tPel','-append');



