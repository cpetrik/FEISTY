% FEISTY output at all locations

clear
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';
harv = 'All_fishobs';

esms = {'IPSL','GFDL','HAD'};

for z=2:3

vers = esms{z};

fpath=['/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/',vers,'/'];

%% SP
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_rec_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[nid,nt] = size(rec);

SP.rec = rec;
SP.gamma = gamma;
clear rec gamma

%% SF
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_rec_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.rec = rec(:,1:nt);
SF.gamma = gamma;
clear rec gamma

% SD
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_rec_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.rec = rec;
SD.gamma = gamma;
clear rec gamma

% MP
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_rec_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.rec = rec;
MP.gamma = gamma;
clear rec gamma

% MF
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_rec_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.rec = rec;
MF.gamma = gamma;
clear rec gamma

% MD
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_rec_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.rec = rec;
MD.gamma = gamma;
clear rec gamma

% LP
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_rec_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.rec = rec;
LP.gamma = gamma;
clear rec gamma

% LD
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_rec_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.rec = rec;
LD.gamma = gamma;
clear rec gamma

%% Take means
%Time
sp_tmrec=mean(SP.rec,1,"omitnan");
sf_tmrec=mean(SF.rec,1,"omitnan");
sd_tmrec=mean(SD.rec,1,"omitnan");
mp_tmrec=mean(MP.rec,1,"omitnan");
mf_tmrec=mean(MF.rec,1,"omitnan");
md_tmrec=mean(MD.rec,1,"omitnan");
lp_tmrec=mean(LP.rec,1,"omitnan");
ld_tmrec=mean(LD.rec,1,"omitnan");

sf_tmgamma=mean(SF.gamma,1,"omitnan");
sp_tmgamma=mean(SP.gamma,1,"omitnan");
sd_tmgamma=mean(SD.gamma,1,"omitnan");
mf_tmgamma=mean(MF.gamma,1,"omitnan");
mp_tmgamma=mean(MP.gamma,1,"omitnan");
md_tmgamma=mean(MD.gamma,1,"omitnan");
lp_tmgamma=mean(LP.gamma,1,"omitnan");
ld_tmgamma=mean(LD.gamma,1,"omitnan");

%% All years
%lyr=time((end-12+1):end);
%lyr=1:12;
sp_smrec=mean(SP.rec,2,"omitnan");
sf_smrec=mean(SF.rec,2,"omitnan");
sd_smrec=mean(SD.rec,2,"omitnan");
mp_smrec=mean(MP.rec,2,"omitnan");
mf_smrec=mean(MF.rec,2,"omitnan");
md_smrec=mean(MD.rec,2,"omitnan");
lp_smrec=mean(LP.rec,2,"omitnan");
ld_smrec=mean(LD.rec,2,"omitnan");

sf_smgamma=mean(SF.gamma,2,"omitnan");
sp_smgamma=mean(SP.gamma,2,"omitnan");
sd_smgamma=mean(SD.gamma,2,"omitnan");
mf_smgamma=mean(MF.gamma,2,"omitnan");
mp_smgamma=mean(MP.gamma,2,"omitnan");
md_smgamma=mean(MD.gamma,2,"omitnan");
lp_smgamma=mean(LP.gamma,2,"omitnan");
ld_smgamma=mean(LD.gamma,2,"omitnan");

%% Each year
a = 1:12:nt; % start of each yr
b = 12:12:nt; % end of each yr
ymrB = NaN*ones(length(sf_smrec),(nt/12));
ymrSF = ymrB;
ymrSP = ymrB;
ymrSD = ymrB;
ymrMF = ymrB;
ymrMP = ymrB;
ymrMD = ymrB;
ymrLP = ymrB;
ymrLD = ymrB;
ymgSF = ymrB;
ymgSP = ymrB;
ymgSD = ymrB;
ymgMF = ymrB;
ymgMP = ymrB;
ymgMD = ymrB;
ymgLP = ymrB;
ymgLD = ymrB;
for i = 1:(nt/12)
    %yr = (i+1987);
    %! Put vars of netcdf file
    ymrSF(:,i) = mean(SF.rec(:,a(i):b(i)),2,"omitnan");
    ymrSP(:,i) = mean(SP.rec(:,a(i):b(i)),2,"omitnan");
    ymrSD(:,i) = mean(SD.rec(:,a(i):b(i)),2,"omitnan");
    ymrMF(:,i) = mean(MF.rec(:,a(i):b(i)),2,"omitnan");
    ymrMP(:,i) = mean(MP.rec(:,a(i):b(i)),2,"omitnan");
    ymrMD(:,i) = mean(MD.rec(:,a(i):b(i)),2,"omitnan");
    ymrLP(:,i) = mean(LP.rec(:,a(i):b(i)),2,"omitnan");
    ymrLD(:,i) = mean(LD.rec(:,a(i):b(i)),2,"omitnan");

    ymgSF(:,i) = mean(SF.gamma(:,a(i):b(i)),2,"omitnan");
    ymgSP(:,i) = mean(SP.gamma(:,a(i):b(i)),2,"omitnan");
    ymgSD(:,i) = mean(SD.gamma(:,a(i):b(i)),2,"omitnan");
    ymgMF(:,i) = mean(MF.gamma(:,a(i):b(i)),2,"omitnan");
    ymgMP(:,i) = mean(MP.gamma(:,a(i):b(i)),2,"omitnan");
    ymgMD(:,i) = mean(MD.gamma(:,a(i):b(i)),2,"omitnan");
    ymgLP(:,i) = mean(LP.gamma(:,a(i):b(i)),2,"omitnan");
    ymgLD(:,i) = mean(LD.gamma(:,a(i):b(i)),2,"omitnan");
end

%%
save([fpath 'Means_Project_' vers '_' harv '_rec_' cfile '.mat'],...
    'sf_smrec','sp_smrec','sd_smrec','mf_smrec','mp_smrec','md_smrec',...
    'lp_smrec','ld_smrec',...
    'sf_tmrec','sp_tmrec','sd_tmrec','mf_tmrec','mp_tmrec','md_tmrec',...
    'lp_tmrec','ld_tmrec',...
    'sf_tmgamma','sp_tmgamma','sd_tmgamma','mf_tmgamma','mp_tmgamma','md_tmgamma',...
    'lp_tmgamma','ld_tmgamma',...
    'sf_smgamma','sp_smgamma','sd_smgamma','mf_smgamma','mp_smgamma','md_smgamma',...
    'lp_smgamma','ld_smgamma',...
    'ymrSF','ymrSP','ymrSD','ymrMF','ymrMP','ymrMD','ymrLP','ymrLD',...
    'ymgSF','ymgSP','ymgSD','ymgMF','ymgMP','ymgMD','ymgLP','ymgLD');


end
