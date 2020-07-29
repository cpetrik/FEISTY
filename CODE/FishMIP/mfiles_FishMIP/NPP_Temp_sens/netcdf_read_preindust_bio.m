% FEISTY output at all locations

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

fpath=['/Volumes/FEISTY/NC/FishMIP/CESM1-BEC/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'Preindust_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(biomass);

SP.bio = biomass;
clear biomass

% SF
ncid = netcdf.open([fpath 'Preindust_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
clear biomass 

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
clear biomass 

% MP
ncid = netcdf.open([fpath 'Preindust_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
clear biomass

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
clear biomass

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
clear biomass

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
clear biomass

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
clear biomass

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
clear biomass

%% Take means and totals
%Total system carbon biomass, tsb, gCm-2, All primary producers and consumers 
%Total consumer carbon biomass density, tcb, gCm-2 All consumers (trophic level >1, vertebrates and invertebrates)
%Carbon biomass density of consumers > 10cm, b10cm, gCm-2 If asymptotic length (Linf) is > 10cm, include in > 10cm class  
%Carbon biomass density of consumers > 30cm, b30cm, gCm-2 If asymptotic length (Linf) is > 30cm, include in > 30cm class
%Monthly or annual 

MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
t=time;
mo=t/12;
mo=mo+1850;
yr18=find(mo>=1851 & mo<1901);  %1851-1900
yr19=find(mo>=1951 & mo<2001);  %1951-2000
yr20=find(mo>=2051 & mo<2101);  %2051-2100

% Ryan comparing 2090-2100 to 1860-1870
yrP=find(mo>=1860 & mo<1870);  
yrF=find(mo>=2090 & mo<2100);  

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

%% Means 
%y = 1850+(1/12):(1/12):2101;
sp_mean18=mean(SP.bio(:,yr18),2);
sf_mean18=mean(SF.bio(:,yr18),2);
sd_mean18=mean(SD.bio(:,yr18),2);
mp_mean18=mean(MP.bio(:,yr18),2);
mf_mean18=mean(MF.bio(:,yr18),2);
md_mean18=mean(MD.bio(:,yr18),2);
lp_mean18=mean(LP.bio(:,yr18),2);
ld_mean18=mean(LD.bio(:,yr18),2);
b_mean18 =mean(Bent.bio(:,yr18),2);

sp_mean19=mean(SP.bio(:,yr19),2);
sf_mean19=mean(SF.bio(:,yr19),2);
sd_mean19=mean(SD.bio(:,yr19),2);
mp_mean19=mean(MP.bio(:,yr19),2);
mf_mean19=mean(MF.bio(:,yr19),2);
md_mean19=mean(MD.bio(:,yr19),2);
lp_mean19=mean(LP.bio(:,yr19),2);
ld_mean19=mean(LD.bio(:,yr19),2);
b_mean19 =mean(Bent.bio(:,yr19),2);

sp_mean20=mean(SP.bio(:,yr20),2);
sf_mean20=mean(SF.bio(:,yr20),2);
sd_mean20=mean(SD.bio(:,yr20),2);
mp_mean20=mean(MP.bio(:,yr20),2);
mf_mean20=mean(MF.bio(:,yr20),2);
md_mean20=mean(MD.bio(:,yr20),2);
lp_mean20=mean(LP.bio(:,yr20),2);
ld_mean20=mean(LD.bio(:,yr20),2);
b_mean20 =mean(Bent.bio(:,yr20),2);

sp_meanP=mean(SP.bio(:,yrP),2);
sf_meanP=mean(SF.bio(:,yrP),2);
sd_meanP=mean(SD.bio(:,yrP),2);
mp_meanP=mean(MP.bio(:,yrP),2);
mf_meanP=mean(MF.bio(:,yrP),2);
md_meanP=mean(MD.bio(:,yrP),2);
lp_meanP=mean(LP.bio(:,yrP),2);
ld_meanP=mean(LD.bio(:,yrP),2);
b_meanP =mean(Bent.bio(:,yrP),2);

sp_meanF=mean(SP.bio(:,yrF),2);
sf_meanF=mean(SF.bio(:,yrF),2);
sd_meanF=mean(SD.bio(:,yrF),2);
mp_meanF=mean(MP.bio(:,yrF),2);
mf_meanF=mean(MF.bio(:,yrF),2);
md_meanF=mean(MD.bio(:,yrF),2);
lp_meanF=mean(LP.bio(:,yrF),2);
ld_meanF=mean(LD.bio(:,yrF),2);
b_meanF =mean(Bent.bio(:,yrF),2);

%% Annual totals
st=1:12:length(time);
en=12:12:length(time);
for m=1:length(en)
    yr1 = st(m):en(m);
    sp_tot(:,m)=sum(SP.bio(:,yr1).*MNTH,2);
    sf_tot(:,m)=sum(SF.bio(:,yr1).*MNTH,2);
    sd_tot(:,m)=sum(SD.bio(:,yr1).*MNTH,2);
    mp_tot(:,m)=sum(MP.bio(:,yr1).*MNTH,2);
    mf_tot(:,m)=sum(MF.bio(:,yr1).*MNTH,2);
    md_tot(:,m)=sum(MD.bio(:,yr1).*MNTH,2);
    lp_tot(:,m)=sum(LP.bio(:,yr1).*MNTH,2);
    ld_tot(:,m)=sum(LD.bio(:,yr1).*MNTH,2);
    b_tot(:,m)=sum(Bent.bio(:,yr1).*MNTH,2);

    sp_mean(:,m)=mean(SP.bio(:,yr1),2);
    sf_mean(:,m)=mean(SF.bio(:,yr1),2);
    sd_mean(:,m)=mean(SD.bio(:,yr1),2);
    mp_mean(:,m)=mean(MP.bio(:,yr1),2);
    mf_mean(:,m)=mean(MF.bio(:,yr1),2);
    md_mean(:,m)=mean(MD.bio(:,yr1),2);
    lp_mean(:,m)=mean(LP.bio(:,yr1),2);
    ld_mean(:,m)=mean(LD.bio(:,yr1),2);
    b_mean(:,m)=mean(Bent.bio(:,yr1),2);
end

%%
save([fpath 'Means_Preindust_' cfile '.mat'],'time','mo',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean',...
    'sf_tot','sp_tot','sd_tot',...
    'mf_tot','mp_tot','md_tot',...
    'lp_tot','ld_tot','b_tot',...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean',...
    'sf_mean18','sp_mean18','sd_mean18',...
    'mf_mean18','mp_mean18','md_mean18',...
    'lp_mean18','ld_mean18','b_mean18',...
    'sf_mean19','sp_mean19','sd_mean19',...
    'mf_mean19','mp_mean19','md_mean19',...
    'lp_mean19','ld_mean19','b_mean19',...
    'sf_mean20','sp_mean20','sd_mean20',...
    'mf_mean20','mp_mean20','md_mean20',...
    'lp_mean20','ld_mean20','b_mean20',...
    'sf_meanP','sp_meanP','sd_meanP',...
    'mf_meanP','mp_meanP','md_meanP',...
    'lp_meanP','ld_meanP','b_meanP',...
    'sf_meanF','sp_meanF','sd_meanF',...
    'mf_meanF','mp_meanF','md_meanF',...
    'lp_meanF','ld_meanF','b_meanF');







