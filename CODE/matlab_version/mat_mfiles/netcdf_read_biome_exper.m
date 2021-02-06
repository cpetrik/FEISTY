% CORE-forced climatology run

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
exper = 'Biome_exper_control_';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/' exper];

%% SP
ncid = netcdf.open([fpath 'All_fish03_sml_p.nc'],'NC_NOWRITE');
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

% SF
ncid = netcdf.open([fpath 'All_fish03_sml_f.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'All_fish03_sml_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'All_fish03_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
Med_p.bio = biomass(:,nt);

clear biomass 

% MF
ncid = netcdf.open([fpath 'All_fish03_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
Med_f.bio = biomass(:,nt);

clear biomass

% MD
ncid = netcdf.open([fpath 'All_fish03_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
Med_d.bio = biomass(:,nt);

clear biomass 

% LP
ncid = netcdf.open([fpath 'All_fish03_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
Lrg_p.bio = biomass(:,nt);

clear biomass 

% LD
ncid = netcdf.open([fpath 'All_fish03_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
Lrg_d.bio = biomass(:,nt);

clear biomass 

% Benthic material
ncid = netcdf.open([fpath 'All_fish03_bent.nc'],'NC_NOWRITE');
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

% take means of each slope
% first 20 rows all slope = 0 in biome 1
% next 20 rows all slope = -0.3 in biome 1
% 221:240 all slope = 0 in biome 2
% etc.
st=1:20:nid;
en=20:20:nid;

for n=1:length(st)
    sp_mean(n,:)=nanmean(SP.bio(st(n):en(n),:),1);
    sf_mean(n,:)=nanmean(SF.bio(st(n):en(n),:),1);
    sd_mean(n,:)=nanmean(SD.bio(st(n):en(n),:),1);
    mp_mean(n,:)=nanmean(MP.bio(st(n):en(n),:),1);
    mf_mean(n,:)=nanmean(MF.bio(st(n):en(n),:),1);
    md_mean(n,:)=nanmean(MD.bio(st(n):en(n),:),1);
    lp_mean(n,:)=nanmean(LP.bio(st(n):en(n),:),1);
    ld_mean(n,:)=nanmean(LD.bio(st(n):en(n),:),1);
    b_mean(n,:)=nanmean(Bent.bio(st(n):en(n),:),1);
    
end

% Means of each biome
st=1:220:nid;
en=220:220:nid;

for n=1:length(st)
    sp_meanB(n,:)=nanmean(SP.bio(st(n):en(n),:),1);
    sf_meanB(n,:)=nanmean(SF.bio(st(n):en(n),:),1);
    sd_meanB(n,:)=nanmean(SD.bio(st(n):en(n),:),1);
    mp_meanB(n,:)=nanmean(MP.bio(st(n):en(n),:),1);
    mf_meanB(n,:)=nanmean(MF.bio(st(n):en(n),:),1);
    md_meanB(n,:)=nanmean(MD.bio(st(n):en(n),:),1);
    lp_meanB(n,:)=nanmean(LP.bio(st(n):en(n),:),1);
    ld_meanB(n,:)=nanmean(LD.bio(st(n):en(n),:),1);
    b_meanB(n,:)=nanmean(Bent.bio(st(n):en(n),:),1);
    
end

%% Save means
save([fpath 'Means_fished_' cfile '.mat'],'time',...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'b_mean','lp_mean','ld_mean',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'b_tmean','lp_tmean','ld_tmean',...
    'sf_meanB','sp_meanB','sd_meanB',...
    'mf_meanB','mp_meanB','md_meanB',...
    'lp_meanB','ld_meanB','b_meanB');
    

