% FEISTY output at all locations

clear
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';
vers = 'IPSL';
harv = 'All_fishobs';
%Hist_IPSL_All_fishobs_bent.nc

fpath=['/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/',vers,'/'];

%% Benthic material
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

nt = length(time);

Bent.bio = biomass;
BENT.bio = biomass(:,nt);
clear biomass

% SP
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.bio = biomass;
Sml_p.bio = biomass(:,nt);
clear biomass

% SF
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_sml_f.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_sml_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
MP.yield = yield;
Med_p.bio = biomass(:,nt);
clear biomass yield

% MF
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
MF.yield = yield;
Med_f.bio = biomass(:,nt);
clear biomass yield

% MD
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
MD.yield = yield;
Med_d.bio = biomass(:,nt);
clear biomass yield

% LP
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
LP.yield = yield;
Lrg_p.bio = biomass(:,nt);
clear biomass yield

% LD
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
LD.yield = yield;
Lrg_d.bio = biomass(:,nt);
clear biomass yield

%% Take means
%Time
sp_tmean=mean(SP.bio,1,"omitnan");
sf_tmean=mean(SF.bio,1,"omitnan");
sd_tmean=mean(SD.bio,1,"omitnan");
mp_tmean=mean(MP.bio,1,"omitnan");
mf_tmean=mean(MF.bio,1,"omitnan");
md_tmean=mean(MD.bio,1,"omitnan");
lp_tmean=mean(LP.bio,1,"omitnan");
ld_tmean=mean(LD.bio,1,"omitnan");
b_tmean=mean(Bent.bio,1,"omitnan");

mf_tmy=mean(MF.yield,1,"omitnan");
mp_tmy=mean(MP.yield,1,"omitnan");
md_tmy=mean(MD.yield,1,"omitnan");
lp_tmy=mean(LP.yield,1,"omitnan");
ld_tmy=mean(LD.yield,1,"omitnan");

%% All years
%lyr=time((end-12+1):end);
%lyr=1:12;
sp_mean=mean(SP.bio,2,"omitnan");
sf_mean=mean(SF.bio,2,"omitnan");
sd_mean=mean(SD.bio,2,"omitnan");
mp_mean=mean(MP.bio,2,"omitnan");
mf_mean=mean(MF.bio,2,"omitnan");
md_mean=mean(MD.bio,2,"omitnan");
lp_mean=mean(LP.bio,2,"omitnan");
ld_mean=mean(LD.bio,2,"omitnan");
b_mean=mean(Bent.bio,2,"omitnan");

mf_my=mean(MF.yield,2,"omitnan");
mp_my=mean(MP.yield,2,"omitnan");
md_my=mean(MD.yield,2,"omitnan");
lp_my=mean(LP.yield,2,"omitnan");
ld_my=mean(LD.yield,2,"omitnan");

%% Each year
a = 1:12:nt; % start of each yr
b = 12:12:nt; % end of each yr
mB = NaN*ones(length(b_mean),(nt/12));
mSF = mB;
mSP = mB;
mSD = mB;
mMF = mB;
mMP = mB;
mMD = mB;
mLP = mB;
mLD = mB;
for i = 1:(nt/12)
    %yr = (i+1987);
    %! Put vars of netcdf file
    mB(:,i) = mean(Bent.bio(:,a(i):b(i)),2,"omitnan");
    mSF(:,i) = mean(SF.bio(:,a(i):b(i)),2,"omitnan");
    mSP(:,i) = mean(SP.bio(:,a(i):b(i)),2,"omitnan");
    mSD(:,i) = mean(SD.bio(:,a(i):b(i)),2,"omitnan");
    mMF(:,i) = mean(MF.bio(:,a(i):b(i)),2,"omitnan");
    mMP(:,i) = mean(MP.bio(:,a(i):b(i)),2,"omitnan");
    mMD(:,i) = mean(MD.bio(:,a(i):b(i)),2,"omitnan");
    mLP(:,i) = mean(LP.bio(:,a(i):b(i)),2,"omitnan");
    mLD(:,i) = mean(LD.bio(:,a(i):b(i)),2,"omitnan");
end

%%
save([fpath 'Means_Hist_' vers '_' harv '_' cfile '.mat'],...
    'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
    'lp_mean','ld_mean','sf_tmean','sp_tmean','sd_tmean','mf_tmean','mp_tmean',...
    'md_tmean','b_tmean','lp_tmean','ld_tmean','time',...
    'mf_tmy','mp_tmy','md_tmy','lp_tmy','ld_tmy',...
    'mf_my','mp_my','md_my','lp_my','ld_my',...
    'mB','mSF','mSP','mSD','mMF','mMP','mMD','mLP','mLD');

% Save last year for initializing forecast runs
save([fpath 'Last_mo_Hist_' vers '_' harv '_' cfile '.mat'],'Sml_f','Sml_p','Sml_d',...
    'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')




