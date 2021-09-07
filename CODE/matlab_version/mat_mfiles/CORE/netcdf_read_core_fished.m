% CORE-forced run

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'Core_All_fish03_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[nid,nt] = size(biomass);
SP.bio = biomass;
Sml_p.bio = biomass(:,nt);

clear biomass 

%% SF
ncid = netcdf.open([fpath 'Core_All_fish03_sml_f.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Core_All_fish03_sml_d.nc'],'NC_NOWRITE');
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

%% MP
ncid = netcdf.open([fpath 'Core_All_fish03_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
Med_p.bio = biomass(:,nt);
MP.yield = yield;

clear biomass yield

% MF
ncid = netcdf.open([fpath 'Core_All_fish03_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
Med_f.bio = biomass(:,nt);
MF.yield = yield;

clear biomassyield

% MD
ncid = netcdf.open([fpath 'Core_All_fish03_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
Med_d.bio = biomass(:,nt);
MD.yield = yield;

clear biomass yield

% LP
ncid = netcdf.open([fpath 'Core_All_fish03_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
Lrg_p.bio = biomass(:,nt);
LP.yield = yield;

clear biomass yield

% LD
ncid = netcdf.open([fpath 'Core_All_fish03_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
Lrg_d.bio = biomass(:,nt);
LD.yield = yield;

clear biomass yield

% Benthic material
ncid = netcdf.open([fpath 'Core_All_fish03_bent.nc'],'NC_NOWRITE');
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
nt = length(time);
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];


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

%% Each year
st=1:12:length(time);
en=12:12:length(time);

for n=1:length(st)
    sp_mean(:,n)=nanmean(SP.bio(:,st(n):en(n)),2);
    sf_mean(:,n)=nanmean(SF.bio(:,st(n):en(n)),2);
    sd_mean(:,n)=nanmean(SD.bio(:,st(n):en(n)),2);
    mp_mean(:,n)=nanmean(MP.bio(:,st(n):en(n)),2);
    mf_mean(:,n)=nanmean(MF.bio(:,st(n):en(n)),2);
    md_mean(:,n)=nanmean(MD.bio(:,st(n):en(n)),2);
    lp_mean(:,n)=nanmean(LP.bio(:,st(n):en(n)),2);
    ld_mean(:,n)=nanmean(LD.bio(:,st(n):en(n)),2);
    b_mean(:,n)=nanmean(Bent.bio(:,st(n):en(n)),2);
    
%     mp_my(:,n)=nanmean(MP.yield(:,st(n):en(n)),2);
%     mf_my(:,n)=nanmean(MF.yield(:,st(n):en(n)),2);
%     md_my(:,n)=nanmean(MD.yield(:,st(n):en(n)),2);
%     lp_my(:,n)=nanmean(LP.yield(:,st(n):en(n)),2);
%     ld_my(:,n)=nanmean(LD.yield(:,st(n):en(n)),2);
end

%% Whole time period mean
sp_mean58=mean(SP.bio,2);
sf_mean58=mean(SF.bio,2);
sd_mean58=mean(SD.bio,2);
mp_mean58=mean(MP.bio,2);
mf_mean58=mean(MF.bio,2);
md_mean58=mean(MD.bio,2);
lp_mean58=mean(LP.bio,2);
ld_mean58=mean(LD.bio,2);
b_mean58=mean(Bent.bio,2);

mf_my58=mean(MF.yield,2);
mp_my58=mean(MP.yield,2);
md_my58=mean(MD.yield,2);
lp_my58=mean(LP.yield,2);
ld_my58=mean(LD.yield,2);

%% Whole time period each month
sp_bio=SP.bio;
sf_bio=SF.bio;
sd_bio=SD.bio;
mp_bio=MP.bio;
mf_bio=MF.bio;
md_bio=MD.bio;
lp_bio=LP.bio;
ld_bio=LD.bio;
b_bio=Bent.bio;

mf_yield=MF.yield;
mp_yield=MP.yield;
md_yield=MD.yield;
lp_yield=LP.yield;
ld_yield=LD.yield;

%% Save means
save([fpath 'Means_core_fished_' cfile '.mat'],'time',...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'b_mean','lp_mean','ld_mean',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'b_tmean','lp_tmean','ld_tmean',...
    'sf_mean58','sp_mean58','sd_mean58',...
    'mf_mean58','mp_mean58','md_mean58',...
    'lp_mean58','ld_mean58','b_mean58',...
    'mf_my58','mp_my58','md_my58',...
    'lp_my58','ld_my58',...
    'sf_bio','sp_bio','sd_bio',...
    'mf_bio','mp_bio','md_bio',...
    'b_bio','lp_bio','ld_bio',...
    'mf_yield','mp_yield','md_yield',...
    'lp_yield','ld_yield');
    

