% POEM output at all locations

clear all
close all

%cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'Climatol_pristine_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.bio = biomass;
SP.prod = prod;
clear biomass prod

% SF
ncid = netcdf.open([fpath 'Climatol_pristine_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass;
SF.prod = prod;
clear biomass prod

% SD
ncid = netcdf.open([fpath 'Climatol_pristine_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
SD.prod = prod;
clear biomass prod

% MP
ncid = netcdf.open([fpath 'Climatol_pristine_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
MP.prod = prod;
clear biomass prod 

%% MF
ncid = netcdf.open([fpath 'Climatol_pristine_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
MF.prod = prod;
clear biomass prod 

% MD
ncid = netcdf.open([fpath 'Climatol_pristine_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
MD.prod = prod;
clear biomass prod 

% LP
ncid = netcdf.open([fpath 'Climatol_pristine_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
LP.prod = prod;
clear biomass prod 

% LD
ncid = netcdf.open([fpath 'Climatol_pristine_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
LD.prod = prod;
clear biomass prod 

% Benthic material
ncid = netcdf.open([fpath 'Climatol_pristine_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

BENT.bio = biomass;
clear biomass

%% Take means and totals
nt = length(time);
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%Time
sp_tmean=mean(SP.bio,1);
sf_tmean=mean(SF.bio(:,1:nt),1);
sd_tmean=mean(SD.bio,1);
mp_tmean=mean(MP.bio,1);
mf_tmean=mean(MF.bio,1);
md_tmean=mean(MD.bio,1);
lp_tmean=mean(LP.bio,1);
ld_tmean=mean(LD.bio,1);
b_tmean=mean(BENT.bio,1);

sp_tmprod=mean(SP.prod,1);
sf_tmprod=mean(SF.prod(:,1:nt),1);
sd_tmprod=mean(SD.prod,1);
mp_tmprod=mean(MP.prod,1);
mf_tmprod=mean(MF.prod,1);
md_tmprod=mean(MD.prod,1);
lp_tmprod=mean(LP.prod,1);
ld_tmprod=mean(LD.prod,1);

% Last year
lyr=time((end-12+1):end);
%Means
sp_mean=mean(SP.bio(:,lyr),2);
sf_mean=mean(SF.bio(:,lyr),2);
sd_mean=mean(SD.bio(:,lyr),2);
mp_mean=mean(MP.bio(:,lyr),2);
mf_mean=mean(MF.bio(:,lyr),2);
md_mean=mean(MD.bio(:,lyr),2);
lp_mean=mean(LP.bio(:,lyr),2);
ld_mean=mean(LD.bio(:,lyr),2);
b_mean=mean(BENT.bio(:,lyr),2);
all_mean1=mean(SP.bio(:,lyr),2)+mean(SF.bio(:,lyr),2)+mean(SD.bio(:,lyr),2)+...
    mean(MP.bio(:,lyr),2)+mean(MF.bio(:,lyr),2)+mean(MD.bio(:,lyr),2)+...
    mean(LP.bio(:,lyr),2)+mean(LD.bio(:,lyr),2);
all_mean2=mean((SP.bio(:,lyr)+SF.bio(:,lyr)+SD.bio(:,lyr)+...
    MP.bio(:,lyr)+MF.bio(:,lyr)+MD.bio(:,lyr)+...
    LP.bio(:,lyr)+LD.bio(:,lyr)),2);

%Medians
all_median1=median(SP.bio(:,lyr),2)+median(SF.bio(:,lyr),2)+median(SD.bio(:,lyr),2)+...
    median(MP.bio(:,lyr),2)+median(MF.bio(:,lyr),2)+median(MD.bio(:,lyr),2)+...
    median(LP.bio(:,lyr),2)+median(LD.bio(:,lyr),2);
all_median2=median((SP.bio(:,lyr)+SF.bio(:,lyr)+SD.bio(:,lyr)+...
    MP.bio(:,lyr)+MF.bio(:,lyr)+MD.bio(:,lyr)+...
    LP.bio(:,lyr)+LD.bio(:,lyr)),2);

sp_mprod=mean(SP.prod(:,lyr),2);
sf_mprod=mean(SF.prod(:,lyr),2);
sd_mprod=mean(SD.prod(:,lyr),2);
mp_mprod=mean(MP.prod(:,lyr),2);
mf_mprod=mean(MF.prod(:,lyr),2);
md_mprod=mean(MD.prod(:,lyr),2);
lp_mprod=mean(LP.prod(:,lyr),2);
ld_mprod=mean(LD.prod(:,lyr),2);

%Totals
sp_tot=sum(SP.bio(:,lyr).*MNTH,2);
sf_tot=sum(SF.bio(:,lyr).*MNTH,2);
sd_tot=sum(SD.bio(:,lyr).*MNTH,2);
mp_tot=sum(MP.bio(:,lyr).*MNTH,2);
mf_tot=sum(MF.bio(:,lyr).*MNTH,2);
md_tot=sum(MD.bio(:,lyr).*MNTH,2);
lp_tot=sum(LP.bio(:,lyr).*MNTH,2);
ld_tot=sum(LD.bio(:,lyr).*MNTH,2);
b_tot=sum(BENT.bio(:,lyr).*MNTH,2);

sp_tprod=sum(SP.prod(:,lyr).*MNTH,2);
sf_tprod=sum(SF.prod(:,lyr).*MNTH,2);
sd_tprod=sum(SD.prod(:,lyr).*MNTH,2);
mp_tprod=sum(MP.prod(:,lyr).*MNTH,2);
mf_tprod=sum(MF.prod(:,lyr).*MNTH,2);
md_tprod=sum(MD.prod(:,lyr).*MNTH,2);
lp_tprod=sum(LP.prod(:,lyr).*MNTH,2);
ld_tprod=sum(LD.prod(:,lyr).*MNTH,2);

%%
save([fpath 'Means_Climatol_pristine_' cfile '.mat'],'time','lyr',...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean',...
    'sf_mprod','sp_mprod','sd_mprod',...
    'mf_mprod','mp_mprod','md_mprod',...
    'lp_mprod','ld_mprod',...
    'sf_tot','sp_tot','sd_tot',...
    'mf_tot','mp_tot','md_tot',...
    'lp_tot','ld_tot','b_tot',...
    'sf_tprod','sp_tprod','sd_tprod',...
    'mf_tprod','mp_tprod','md_tprod',...
    'lp_tprod','ld_tprod',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean',...
    'sf_tmprod','sp_tmprod','sd_tmprod',...
    'mf_tmprod','mp_tmprod','md_tmprod',...
    'lp_tmprod','ld_tmprod',...
    'all_median1','all_median2','all_mean1','all_mean2');%,'-append');

save([fpath 'Means_bio_prod_Climatol_pristine_' cfile '.mat'],...
    'time','lyr',...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean',...
    'sf_mprod','sp_mprod','sd_mprod',...
    'mf_mprod','mp_mprod','md_mprod',...
    'lp_mprod','ld_mprod',...
    'sf_tot','sp_tot','sd_tot',...
    'mf_tot','mp_tot','md_tot',...
    'lp_tot','ld_tot','b_tot',...
    'sf_tprod','sp_tprod','sd_tprod',...
    'mf_tprod','mp_tprod','md_tprod',...
    'lp_tprod','ld_tprod',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean',...
    'sf_tmprod','sp_tmprod','sd_tmprod',...
    'mf_tmprod','mp_tmprod','md_tmprod',...
    'lp_tmprod','ld_tmprod',...
    'all_median1','all_median2','all_mean1','all_mean2');





