% POEM output at all locations

clear all
close all

%cfile = 'Dc_enc70_cmax-metab20_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00400';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];

nt=12*100;

%% SP
ncid = netcdf.open([fpath 'Preindust_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.bio = biomass;
% SP.clev = clev;
% SP.con = con;
% SP.die = die;
% SP.gamma = gamma;
% SP.nu = nu;
SP.prod = prod;
SP.rec = rec;

Sml_p.bio = biomass(:,nt);

clear biomass clev con die gamma nu rec time prod

%% SF
ncid = netcdf.open([fpath 'Preindust_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
% SF.clev = clev(:,1:nt);
% SF.con = con(:,1:nt);
% SF.die = die(:,1:nt);
% SF.gamma = gamma(:,1:nt);
% SF.nu = nu(:,1:nt);
SF.prod = prod(:,1:nt);
SF.rec = rec(:,1:nt);

Sml_f.bio = biomass(:,nt);

clear biomass clev con die gamma nu rec time prod

% SD
ncid = netcdf.open([fpath 'Preindust_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
% SD.clev = clev;
% SD.con = con;
% SD.die = die;
% SD.gamma = gamma;
% SD.nu = nu;
SD.prod = prod;
SD.rec = rec;

Sml_d.bio = biomass(:,nt);

clear biomass clev con die gamma nu rec time prod

%% MP
ncid = netcdf.open([fpath 'Preindust_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
% MP.clev = clev;
% MP.con = con;
% MP.die = die;
% MP.gamma = gamma;
% MP.nu = nu;
MP.prod = prod;
MP.rec = rec;

Med_p.bio = biomass(:,nt);

clear biomass clev con die gamma nu rec time prod

% MF
ncid = netcdf.open([fpath 'Preindust_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
% MF.clev = clev;
% MF.con = con;
% MF.die = die;
% MF.gamma = gamma;
% MF.nu = nu;
MF.prod = prod;
MF.rec = rec;
% MF.rep = rep;

Med_f.bio = biomass(:,nt);

clear biomass clev con die gamma nu rec rep time prod

% MD
ncid = netcdf.open([fpath 'Preindust_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
% MD.clev = clev;
% MD.con = con;
% MD.die = die;
% MD.gamma = gamma;
% MD.nu = nu;
MD.prod = prod;
MD.rec = rec;

Med_d.bio = biomass(:,nt);

clear biomass clev con die gamma nu rec time prod

% LP
ncid = netcdf.open([fpath 'Preindust_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
% LP.clev = clev;
% LP.con = con;
% LP.die = die;
% LP.gamma = gamma;
% LP.nu = nu;
LP.prod = prod;
LP.rec = rec;
% LP.rep = rep;

Lrg_p.bio = biomass(:,nt);

clear biomass clev con die gamma nu rec rep time prod

% LD
ncid = netcdf.open([fpath 'Preindust_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
% LD.clev = clev;
% LD.con = con;
% LD.die = die;
% LD.gamma = gamma;
% LD.nu = nu;
LD.prod = prod;
LD.rec = rec;
% LD.rep = rep;

Lrg_d.bio = biomass(:,nt);

clear biomass clev con die gamma nu rec rep time prod

% Benthic material
ncid = netcdf.open([fpath 'Preindust_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
BENT.bio = biomass(:,nt);
clear biomass 

%% Take means

%Time
sp_tmean=mean(SP.bio,1);
sf_tmean=mean(SF.bio,1);
sd_tmean=mean(SD.bio,1);
mp_tmean=mean(MP.bio,1);
mf_tmean=mean(MF.bio,1);
md_tmean=mean(MD.bio,1);
lp_tmean=mean(LP.bio,1);
ld_tmean=mean(LD.bio,1);
b_tmean=mean(Bent.bio,1);

sp_tprod=mean(SP.prod,1);
sf_tprod=mean(SF.prod,1);
sd_tprod=mean(SD.prod,1);
mp_tprod=mean(MP.prod,1);
mf_tprod=mean(MF.prod,1);
md_tprod=mean(MD.prod,1);
lp_tprod=mean(LP.prod,1);
ld_tprod=mean(LD.prod,1);

sp_trec=mean(SP.rec,1);
sf_trec=mean(SF.rec,1);
sd_trec=mean(SD.rec,1);
mp_trec=mean(MP.rec,1);
mf_trec=mean(MF.rec,1);
md_trec=mean(MD.rec,1);
lp_trec=mean(LP.rec,1);
ld_trec=mean(LD.rec,1);

%% Last 50 years
yr50=time((end-(50*12)+1):end);
sp_mean=mean(SP.bio(:,yr50),2);
sf_mean=mean(SF.bio(:,yr50),2);
sd_mean=mean(SD.bio(:,yr50),2);
mp_mean=mean(MP.bio(:,yr50),2);
mf_mean=mean(MF.bio(:,yr50),2);
md_mean=mean(MD.bio(:,yr50),2);
lp_mean=mean(LP.bio(:,yr50),2);
ld_mean=mean(LD.bio(:,yr50),2);
b_mean=mean(Bent.bio(:,yr50),2);

sp_prod=nanmean(SP.prod(:,yr50),2);
sf_prod=nanmean(SF.prod(:,yr50),2);
sd_prod=nanmean(SD.prod(:,yr50),2);
mp_prod=nanmean(MP.prod(:,yr50),2);
mf_prod=nanmean(MF.prod(:,yr50),2);
md_prod=nanmean(MD.prod(:,yr50),2);
lp_prod=nanmean(LP.prod(:,yr50),2);
ld_prod=nanmean(LD.prod(:,yr50),2);

sp_rec=nanmean(SP.rec(:,yr50),2);
sf_rec=nanmean(SF.rec(:,yr50),2);
sd_rec=nanmean(SD.rec(:,yr50),2);
mp_rec=nanmean(MP.rec(:,yr50),2);
mf_rec=nanmean(MF.rec(:,yr50),2);
md_rec=nanmean(MD.rec(:,yr50),2);
lp_rec=nanmean(LP.rec(:,yr50),2);
ld_rec=nanmean(LD.rec(:,yr50),2);

%% Save means
save([fpath 'Means_preindust_' cfile '.mat'],...
    'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
    'lp_mean','ld_mean','sf_tmean','sp_tmean','sd_tmean','mf_tmean','mp_tmean',...
    'md_tmean','b_tmean','lp_tmean','ld_tmean','time','yr50',...
    'sf_prod','sp_prod','sd_prod','mf_prod','mp_prod',...
    'md_prod','lp_prod','ld_prod','sf_rec','sp_rec','sd_rec','mf_rec',...
    'mp_rec','md_rec','lp_rec','ld_rec',...
    'sf_tprod','sp_tprod','sd_tprod','mf_tprod','mp_tprod',...
    'md_tprod','lp_tprod','ld_tprod','sf_trec','sp_trec','sd_trec','mf_trec',...
    'mp_trec','md_trec','lp_trec','ld_trec');

% Save last year for initializing historical runs
save([fpath 'Last_mo_preindust_' cfile '.mat'],'Sml_f','Sml_p','Sml_d',... 
    'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')


