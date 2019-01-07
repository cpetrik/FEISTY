% POEM output at all locations

clear all
close all


cfile = 'Dc_enc70-b200_m4-b175-k08_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
var = '_con_types';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'Climatol_' harv var '_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.conZ = conZ;
clear conZ 

% SF
ncid = netcdf.open([fpath 'Climatol_' harv var '_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.conZ = conZ;
clear conZ 

% SD
ncid = netcdf.open([fpath 'Climatol_' harv var '_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.conZ = conZ;
clear conZ 

% MP
ncid = netcdf.open([fpath 'Climatol_' harv var '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.conZm = conZm;
MP.conZl = conZl;
MP.conF = conF;
MP.conP = conP;
clear conZm conZl conF conP

%% MF
ncid = netcdf.open([fpath 'Climatol_' harv var '_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.conZm = conZm;
MF.conZl = conZl;
MF.conF = conF;
MF.conP = conP;
clear conZm conZl conF conP

% MD
ncid = netcdf.open([fpath 'Climatol_' harv var '_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.conB = conB;
clear conB

% LP
ncid = netcdf.open([fpath 'Climatol_' harv var '_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.conF = conF;
LP.conP = conP;
clear conF conP

%% LD
ncid = netcdf.open([fpath 'Climatol_' harv var '_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.conF = conF;
LD.conP = conP;
LD.conD = conD;
LD.conB = conB;
clear conB conD conF conP


%% Take means
[ids,nt] = size(LD.conB);

%Time
sp_tmconZm=mean(SP.conZ,1);
sf_tmconZm=mean(SF.conZ,1);
sd_tmconZm=mean(SD.conZ,1);

mp_tmconZm=mean(MP.conZm,1);
mp_tmconZl=mean(MP.conZl,1);
mp_tmconF=mean(MP.conF,1);
mp_tmconP=mean(MP.conP,1);

mf_tmconZm=mean(MF.conZm,1);
mf_tmconZl=mean(MF.conZl,1);
mf_tmconP=mean(MF.conP,1);
mf_tmconF=mean(MF.conF,1);

md_tmconB=mean(MD.conB,1);

lp_tmconF=mean(LP.conF,1);
lp_tmconP=mean(LP.conP,1);

ld_tmconF=mean(LD.conF,1);
ld_tmconP=mean(LD.conP,1);
ld_tmconD=mean(LD.conD,1);
ld_tmconB=mean(LD.conB,1);

%% Last year
time=1:nt;
lyr=time((end-12+1):end);

sp_mconZm=mean(SP.conZ(:,lyr),2);
sf_mconZm=mean(SF.conZ(:,lyr),2);
sd_mconZm=mean(SD.conZ(:,lyr),2);

mp_mconZm=mean(MP.conZm(:,lyr),2);
mp_mconZl=mean(MP.conZl(:,lyr),2);
mp_mconF=mean(MP.conF(:,lyr),2);
mp_mconP=mean(MP.conP(:,lyr),2);

mf_mconZm=mean(MF.conZm(:,lyr),2);
mf_mconZl=mean(MF.conZl(:,lyr),2);
mf_mconF=mean(MF.conF(:,lyr),2);
mf_mconP=mean(MF.conP(:,lyr),2);

md_mconB=mean(MD.conB(:,lyr),2);

lp_mconF=mean(LP.conF(:,lyr),2);
lp_mconP=mean(LP.conP(:,lyr),2);

ld_mconF=mean(LD.conF(:,lyr),2);
ld_mconP=mean(LD.conP(:,lyr),2);
ld_mconD=mean(LD.conD(:,lyr),2);
ld_mconB=mean(LD.conB(:,lyr),2);

%%
save([fpath 'Means_Climatol_' harv '_' cfile '.mat'],'time','lyr',...
    'sp_tmconZm','sf_tmconZm','sd_tmconZm',...
    'mf_tmconZm','mf_tmconZl','mf_tmconF','mf_tmconP',...
    'mp_tmconZm','mp_tmconZl','mp_tmconF','mp_tmconP',...
    'md_tmconB',...
    'lp_tmconF','lp_tmconP',...
    'ld_tmconF','ld_tmconP','ld_tmconD','ld_tmconB',...
    'sp_mconZm','sf_mconZm','sd_mconZm',...
    'mf_mconZm','mf_mconZl','mf_mconF','mf_mconP',...
    'mp_mconZm','mp_mconZl','mp_mconF','mp_mconP',...
    'md_mconB',...
    'lp_mconF','lp_mconP',...
    'ld_mconF','ld_mconP','ld_mconD','ld_mconB',...
    '-append');

save([fpath 'Means',var,'_Climatol_' harv '_' cfile '.mat'],'time','lyr',...
    'sp_tmconZm','sf_tmconZm','sd_tmconZm',...
    'mf_tmconZm','mf_tmconZl','mf_tmconF','mf_tmconP',...
    'mp_tmconZm','mp_tmconZl','mp_tmconF','mp_tmconP',...
    'md_tmconB',...
    'lp_tmconF','lp_tmconP',...
    'ld_tmconF','ld_tmconP','ld_tmconD','ld_tmconB',...
    'sp_mconZm','sf_mconZm','sd_mconZm',...
    'mf_mconZm','mf_mconZl','mf_mconF','mf_mconP',...
    'mp_mconZm','mp_mconZl','mp_mconF','mp_mconP',...
    'md_mconB',...
    'lp_mconF','lp_mconP',...
    'ld_mconF','ld_mconP','ld_mconD','ld_mconB');







