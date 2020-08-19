% FEISTY output at all locations

clear all
close all

cfile = 'Dc_Lam579_enc70-b200_m440-b175-k086_c20-b250_D080_A050_nmort1_BE10_CC80_RE00100';

fpath=['/Volumes/FEISTY/NC/FishMIP/IPSL_CMIP6/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'PreIndust_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%%
[ni,nt] = size(biomass);

SP.bio = biomass;
clear biomass

%% SF
ncid = netcdf.open([fpath 'PreIndust_sml_f.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'PreIndust_sml_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'PreIndust_med_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'PreIndust_med_f.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'PreIndust_med_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'PreIndust_lrg_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'PreIndust_lrg_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'PreIndust_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
clear biomass

%% Take means for my own visualization

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

%% Space
t=time;
mo=t/12;
mo=mo+1850;

yr1=find(mo>1890 & mo<=1900); 
yr2=find(mo>2000 & mo<=2010); 
yr3=find(mo>2090 & mo<=2100); 

sp_mean1=mean(SP.bio(:,yr1),2);
sf_mean1=mean(SF.bio(:,yr1),2);
sd_mean1=mean(SD.bio(:,yr1),2);
mp_mean1=mean(MP.bio(:,yr1),2);
mf_mean1=mean(MF.bio(:,yr1),2);
md_mean1=mean(MD.bio(:,yr1),2);
lp_mean1=mean(LP.bio(:,yr1),2);
ld_mean1=mean(LD.bio(:,yr1),2);
b_mean1 =mean(Bent.bio(:,yr1),2);

sp_mean2=mean(SP.bio(:,yr2),2);
sf_mean2=mean(SF.bio(:,yr2),2);
sd_mean2=mean(SD.bio(:,yr2),2);
mp_mean2=mean(MP.bio(:,yr2),2);
mf_mean2=mean(MF.bio(:,yr2),2);
md_mean2=mean(MD.bio(:,yr2),2);
lp_mean2=mean(LP.bio(:,yr2),2);
ld_mean2=mean(LD.bio(:,yr2),2);
b_mean2 =mean(Bent.bio(:,yr2),2);

sp_mean3=mean(SP.bio(:,yr3),2);
sf_mean3=mean(SF.bio(:,yr3),2);
sd_mean3=mean(SD.bio(:,yr3),2);
mp_mean3=mean(MP.bio(:,yr3),2);
mf_mean3=mean(MF.bio(:,yr3),2);
md_mean3=mean(MD.bio(:,yr3),2);
lp_mean3=mean(LP.bio(:,yr3),2);
ld_mean3=mean(LD.bio(:,yr3),2);
b_mean3 =mean(Bent.bio(:,yr3),2);

save([fpath 'Means_PreIndust_' cfile '.mat'],'time','mo',...
    'yr1','yr2','yr3',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean',...
    'sf_mean1','sp_mean1','sd_mean1',...
    'mf_mean1','mp_mean1','md_mean1',...
    'lp_mean1','ld_mean1','b_mean1',...
    'sf_mean2','sp_mean2','sd_mean2',...
    'mf_mean2','mp_mean2','md_mean2',...
    'lp_mean2','ld_mean2','b_mean2',...
    'sf_mean3','sp_mean3','sd_mean3',...
    'mf_mean3','mp_mean3','md_mean3',...
    'lp_mean3','ld_mean3','b_mean3')

figure
plot(mo,log10(lp_tmean),'b'); hold on;
plot(mo,log10(mf_tmean),'r'); hold on;
plot(mo,log10(ld_tmean),'k'); hold on;

%% Fish-MIP OUTPUTS =================================================
