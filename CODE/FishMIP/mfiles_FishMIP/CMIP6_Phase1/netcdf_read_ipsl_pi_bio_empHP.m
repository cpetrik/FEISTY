% FEISTY output at all locations

clear all
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/MIP/NC/FishMIP/IPSL_CMIP6/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'PreIndust_empHP_sml_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'PreIndust_empHP_sml_f.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'PreIndust_empHP_sml_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'PreIndust_empHP_med_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'PreIndust_empHP_med_f.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'PreIndust_empHP_med_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'PreIndust_empHP_lrg_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'PreIndust_empHP_lrg_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'PreIndust_empHP_bent.nc'],'NC_NOWRITE');
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
sp_tmean=nanmean(SP.bio,1);
sf_tmean=nanmean(SF.bio,1);
sd_tmean=nanmean(SD.bio,1);
mp_tmean=nanmean(MP.bio,1);
mf_tmean=nanmean(MF.bio,1);
md_tmean=nanmean(MD.bio,1);
lp_tmean=nanmean(LP.bio,1);
ld_tmean=nanmean(LD.bio,1);
b_tmean=nanmean(Bent.bio,1);

%% Space
t=time;
mo=t/12;
mo=mo+1950;

yr1=find(mo>1890 & mo<=1900); 
yr2=find(mo>2000 & mo<=2010); 
yr3=find(mo>2090 & mo<=2100); 

sp_mean1=nanmean(SP.bio(:,yr1),2);
sf_mean1=nanmean(SF.bio(:,yr1),2);
sd_mean1=nanmean(SD.bio(:,yr1),2);
mp_mean1=nanmean(MP.bio(:,yr1),2);
mf_mean1=nanmean(MF.bio(:,yr1),2);
md_mean1=nanmean(MD.bio(:,yr1),2);
lp_mean1=nanmean(LP.bio(:,yr1),2);
ld_mean1=nanmean(LD.bio(:,yr1),2);
b_mean1 =nanmean(Bent.bio(:,yr1),2);

sp_mean2=nanmean(SP.bio(:,yr2),2);
sf_mean2=nanmean(SF.bio(:,yr2),2);
sd_mean2=nanmean(SD.bio(:,yr2),2);
mp_mean2=nanmean(MP.bio(:,yr2),2);
mf_mean2=nanmean(MF.bio(:,yr2),2);
md_mean2=nanmean(MD.bio(:,yr2),2);
lp_mean2=nanmean(LP.bio(:,yr2),2);
ld_mean2=nanmean(LD.bio(:,yr2),2);
b_mean2 =nanmean(Bent.bio(:,yr2),2);

sp_mean3=nanmean(SP.bio(:,yr3),2);
sf_mean3=nanmean(SF.bio(:,yr3),2);
sd_mean3=nanmean(SD.bio(:,yr3),2);
mp_mean3=nanmean(MP.bio(:,yr3),2);
mf_mean3=nanmean(MF.bio(:,yr3),2);
md_mean3=nanmean(MD.bio(:,yr3),2);
lp_mean3=nanmean(LP.bio(:,yr3),2);
ld_mean3=nanmean(LD.bio(:,yr3),2);
b_mean3 =nanmean(Bent.bio(:,yr3),2);

save([fpath 'Means_PreIndust_empHP_' cfile '.mat'],'time','mo',...
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
%%
figure
plot(mo,log10(lp_tmean),'b'); hold on;
plot(mo,log10(mf_tmean),'r'); hold on;
plot(mo,log10(ld_tmean),'k'); hold on;

%% Fish-MIP OUTPUTS =================================================

% PREFERRED (all units = gWW/m2)

%total pelagic biomass tpb = 360x180xMOs
allF = SF.bio + MF.bio;
allP = SP.bio + MP.bio + LP.bio;
allPel = allF + allP;

%total demersal biomass tdb = 360x180xMOs
allD = SD.bio + MD.bio + LD.bio;

%total consumber biomass tcb = 360x180xMOs
allC = allF + allP + allD + Bent.bio;

%total consumber biomass in log10 bins tcblog10 = 360x180xMOsx6
%(1g, 10g, 100g, 1kg, 10kg, 100kg)
%Small <1g      (0.001-0.5g; 0.02g mean)    (0.46-3.68cm; mean 1.3cm)
%Med 1-100g     (0.5-250g; 11.2g mean)      (3.68-29.24cm; mean 10.4cm)
%Lrg 1-100kg    (250-125000g; 5600g mean)   (29.24-232.08cm; mean 82.4cm)

% SECONDARY

%total pelagic (Linf <30cm) biomass bp30cm = 360x180xMOs
SPel = allF;

%total pelagic (>=30 cm and <90cm) biomass bp30to90cm = 360x180xMOs
%none

%total pelagic (>=90cm) biomass bp90cm = 360x180xMOs
LPel = allP;

%total demersal (Linf <30cm) biomass bd30cm = 360x180xMOs
%none

%total demersal (>=30 cm and <90cm) biomass bd30to90cm = 360x180xMOs
%none

%total demersal (>=90cm) biomass bd90cm = 360x180xMOs
%LDem = allD;

save([fpath 'PreIndust_empHP_fishMIP_outputs_monthly_' cfile '.mat'],'time','mo',...
    'allPel','allD','allC','SPel','LPel');
