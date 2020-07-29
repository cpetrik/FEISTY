% FEISTY output at all locations

clear all
close all

cfile = 'NoDc_enc70-b200_m4-b175-k086_c20-b250_noD_J100_A050_Sm025_nmort1_BE00_noCC_RE00100';

fpath=['/Volumes/FEISTY/NC/FishMIP/CESM1-BEC/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'Historic_pristine_sml_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Historic_pristine_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
clear biomass 

% MP
ncid = netcdf.open([fpath 'Historic_pristine_med_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Historic_pristine_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
clear biomass

% LP
ncid = netcdf.open([fpath 'Historic_pristine_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
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
mp_tmean=mean(MP.bio,1);
mf_tmean=mean(MF.bio,1);
lp_tmean=mean(LP.bio,1);

%% Means 
%y = 1850+(1/12):(1/12):2101;
sp_mean18=mean(SP.bio(:,yr18),2);
sf_mean18=mean(SF.bio(:,yr18),2);
mp_mean18=mean(MP.bio(:,yr18),2);
mf_mean18=mean(MF.bio(:,yr18),2);
lp_mean18=mean(LP.bio(:,yr18),2);

sp_mean19=mean(SP.bio(:,yr19),2);
sf_mean19=mean(SF.bio(:,yr19),2);
mp_mean19=mean(MP.bio(:,yr19),2);
mf_mean19=mean(MF.bio(:,yr19),2);
lp_mean19=mean(LP.bio(:,yr19),2);

sp_mean20=mean(SP.bio(:,yr20),2);
sf_mean20=mean(SF.bio(:,yr20),2);
mp_mean20=mean(MP.bio(:,yr20),2);
mf_mean20=mean(MF.bio(:,yr20),2);
lp_mean20=mean(LP.bio(:,yr20),2);

sp_meanP=mean(SP.bio(:,yrP),2);
sf_meanP=mean(SF.bio(:,yrP),2);
mp_meanP=mean(MP.bio(:,yrP),2);
mf_meanP=mean(MF.bio(:,yrP),2);
lp_meanP=mean(LP.bio(:,yrP),2);

sp_meanF=mean(SP.bio(:,yrF),2);
sf_meanF=mean(SF.bio(:,yrF),2);
mp_meanF=mean(MP.bio(:,yrF),2);
mf_meanF=mean(MF.bio(:,yrF),2);
lp_meanF=mean(LP.bio(:,yrF),2);

%% Annual totals
st=1:12:length(time);
en=12:12:length(time);
for m=1:length(en)
    yr1 = st(m):en(m);
    sp_tot(:,m)=sum(SP.bio(:,yr1).*MNTH,2);
    sf_tot(:,m)=sum(SF.bio(:,yr1).*MNTH,2);
    mp_tot(:,m)=sum(MP.bio(:,yr1).*MNTH,2);
    mf_tot(:,m)=sum(MF.bio(:,yr1).*MNTH,2);
    lp_tot(:,m)=sum(LP.bio(:,yr1).*MNTH,2);
    
    sp_mean(:,m)=mean(SP.bio(:,yr1),2);
    sf_mean(:,m)=mean(SF.bio(:,yr1),2);
    mp_mean(:,m)=mean(MP.bio(:,yr1),2);
    mf_mean(:,m)=mean(MF.bio(:,yr1),2);
    lp_mean(:,m)=mean(LP.bio(:,yr1),2);
    
end

%%
save([fpath 'Means_Historic_pristine_' cfile '.mat'],'time','mo',...
    'sf_tmean','sp_tmean',...
    'mf_tmean','mp_tmean',...
    'lp_tmean',...
    'sf_tot','sp_tot',...
    'mf_tot','mp_tot',...
    'lp_tot',...
    'sf_mean','sp_mean',...
    'mf_mean','mp_mean',...
    'lp_mean',...
    'sf_mean18','sp_mean18',...
    'mf_mean18','mp_mean18',...
    'lp_mean18',...
    'sf_mean19','sp_mean19',...
    'mf_mean19','mp_mean19',...
    'lp_mean19',...
    'sf_mean20','sp_mean20',...
    'mf_mean20','mp_mean20',...
    'lp_mean20',...
    'sf_meanP','sp_meanP',...
    'mf_meanP','mp_meanP',...
    'lp_meanP',...
    'sf_meanF','sp_meanF',...
    'mf_meanF','mp_meanF',...
    'lp_meanF');

%% Save last year for initializing forecast runs
Sml_f.bio = nanmean(SF.bio(:,nt-11:nt),2);
Sml_p.bio = nanmean(SP.bio(:,nt-11:nt),2);
Med_f.bio = nanmean(MF.bio(:,nt-11:nt),2);
Med_p.bio = nanmean(MP.bio(:,nt-11:nt),2);
Lrg_p.bio = nanmean(LP.bio(:,nt-11:nt),2);

save([fpath 'Last_mo_historic_pristine_' cfile '.mat'],'Sml_f','Sml_p',... 
    'Med_f','Med_p','Lrg_p')





