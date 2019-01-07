% FEISTY output CC
% 3 km model

clear all
close all


%cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
%harv = 'pristine';
harv = 'All_fish01';

ad = 10:10:100;
for n=1:length(ad)
J = ad(n);
if J<100
    cfile = ['Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J0',num2str(J),'_A100',...
    '_Sm075_nmort1_BE08_noCC_RE00100'];
else
    cfile = ['Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J',num2str(J),'_A100',...
    '_Sm075_nmort1_BE08_noCC_RE00100'];
end

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/CalCurr/'];

%% Benthic material
ncid = netcdf.open([fpath 'Historic3km_' harv '_bent.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Historic3km_' harv '_sml_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Historic3km_' harv '_sml_f.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Historic3km_' harv '_sml_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Historic3km_' harv '_med_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Historic3km_' harv '_med_f.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Historic3km_' harv '_med_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Historic3km_' harv '_lrg_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Historic3km_' harv '_lrg_d.nc'],'NC_NOWRITE');
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
sp_tmean=mean(SP.bio,1);
sf_tmean=mean(SF.bio(:,1:nt),1);
sd_tmean=mean(SD.bio,1);
mp_tmean=mean(MP.bio,1);
mf_tmean=mean(MF.bio,1);
md_tmean=mean(MD.bio,1);
lp_tmean=mean(LP.bio,1);
ld_tmean=mean(LD.bio,1);
b_tmean=mean(Bent.bio,1);

mf_tmy=mean(MF.yield,1);
mp_tmy=mean(MP.yield,1);
md_tmy=mean(MD.yield,1);
lp_tmy=mean(LP.yield,1);
ld_tmy=mean(LD.yield,1);

%% All years
%lyr=time((end-12+1):end);
%lyr=1:12;
sp_mean=mean(SP.bio,2);
sf_mean=mean(SF.bio,2);
sd_mean=mean(SD.bio,2);
mp_mean=mean(MP.bio,2);
mf_mean=mean(MF.bio,2);
md_mean=mean(MD.bio,2);
lp_mean=mean(LP.bio,2);
ld_mean=mean(LD.bio,2);
b_mean=mean(Bent.bio,2);

mf_my=mean(MF.yield,2);
mp_my=mean(MP.yield,2);
md_my=mean(MD.yield,2);
lp_my=mean(LP.yield,2);
ld_my=mean(LD.yield,2);

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
    mB(:,i) = mean(Bent.bio(:,a(i):b(i)),2);
    mSF(:,i) = mean(SF.bio(:,a(i):b(i)),2);
    mSP(:,i) = mean(SP.bio(:,a(i):b(i)),2);
    mSD(:,i) = mean(SD.bio(:,a(i):b(i)),2);
    mMF(:,i) = mean(MF.bio(:,a(i):b(i)),2);
    mMP(:,i) = mean(MP.bio(:,a(i):b(i)),2);
    mMD(:,i) = mean(MD.bio(:,a(i):b(i)),2);
    mLP(:,i) = mean(LP.bio(:,a(i):b(i)),2);
    mLD(:,i) = mean(LD.bio(:,a(i):b(i)),2);
end

%% Each season
a = 1:3:nt; % start of each 3-mo period
b = 3:3:nt; % end of each 3-mo period
sB = NaN*ones(length(b_mean),(nt/12));
sSF = sB;
sSP = sB;
sSD = sB;
sMF = sB;
sMP = sB;
sMD = sB;
sLP = sB;
sLD = sB;
for i = 1:(nt/3)
    sB(:,i) = mean(Bent.bio(:,a(i):b(i)),2);
    sSF(:,i) = mean(SF.bio(:,a(i):b(i)),2);
    sSP(:,i) = mean(SP.bio(:,a(i):b(i)),2);
    sSD(:,i) = mean(SD.bio(:,a(i):b(i)),2);
    sMF(:,i) = mean(MF.bio(:,a(i):b(i)),2);
    sMP(:,i) = mean(MP.bio(:,a(i):b(i)),2);
    sMD(:,i) = mean(MD.bio(:,a(i):b(i)),2);
    sLP(:,i) = mean(LP.bio(:,a(i):b(i)),2);
    sLD(:,i) = mean(LD.bio(:,a(i):b(i)),2);
end
    
%%
save([fpath 'Means_Historic3km_' harv '_' cfile '.mat'],...
    'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
    'lp_mean','ld_mean','sf_tmean','sp_tmean','sd_tmean','mf_tmean','mp_tmean',...
    'md_tmean','b_tmean','lp_tmean','ld_tmean','time',...
    'mf_tmy','mp_tmy','md_tmy','lp_tmy','ld_tmy',...
    'mf_my','mp_my','md_my','lp_my','ld_my',...
    'mB','mSF','mSP','mSD','mMF','mMP','mMD','mLP','mLD',...
    'sB','sSF','sSP','sSD','sMF','sMP','sMD','sLP','sLD');

% Save last year for initializing forecast runs
% save([fpath 'Last_mo_Historic3km_' harv '_' cfile '.mat'],'Sml_f','Sml_p','Sml_d',... 
%     'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')

end





