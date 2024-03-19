% FEISTY output at all locations

clear
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';
vers = 'HAD';
harv = 'All_fishobs';
%Hist_IPSL_All_fishobs_bent.nc

fpath=['/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/',vers,'/'];

%% SP
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_prod_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[nid,nt] = size(prod);

SP.prod = prod;
SP.nu = nu;
clear prod nu

%% SF
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_prod_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.prod = prod(:,1:nt);
SF.nu = nu;
clear prod nu

% SD
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_prod_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.prod = prod;
SD.nu = nu;
clear prod nu

% MP
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_prod_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.prod = prod;
MP.nu = nu;
clear prod nu

% MF
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_prod_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.prod = prod;
MF.nu = nu;
clear prod nu

% MD
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_prod_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.prod = prod;
MD.nu = nu;
clear prod nu

% LP
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_prod_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.prod = prod;
LP.nu = nu;
clear prod nu

% LD
ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_prod_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.prod = prod;
LD.nu = nu;
clear prod nu

%% Take means
%Time
sp_tmprod=mean(SP.prod,1,"omitnan");
sf_tmprod=mean(SF.prod,1,"omitnan");
sd_tmprod=mean(SD.prod,1,"omitnan");
mp_tmprod=mean(MP.prod,1,"omitnan");
mf_tmprod=mean(MF.prod,1,"omitnan");
md_tmprod=mean(MD.prod,1,"omitnan");
lp_tmprod=mean(LP.prod,1,"omitnan");
ld_tmprod=mean(LD.prod,1,"omitnan");

sf_tmnu=mean(SF.nu,1,"omitnan");
sp_tmnu=mean(SP.nu,1,"omitnan");
sd_tmnu=mean(SD.nu,1,"omitnan");
mf_tmnu=mean(MF.nu,1,"omitnan");
mp_tmnu=mean(MP.nu,1,"omitnan");
md_tmnu=mean(MD.nu,1,"omitnan");
lp_tmnu=mean(LP.nu,1,"omitnan");
ld_tmnu=mean(LD.nu,1,"omitnan");

%% All years
%lyr=time((end-12+1):end);
%lyr=1:12;
sp_smprod=mean(SP.prod,2,"omitnan");
sf_smprod=mean(SF.prod,2,"omitnan");
sd_smprod=mean(SD.prod,2,"omitnan");
mp_smprod=mean(MP.prod,2,"omitnan");
mf_smprod=mean(MF.prod,2,"omitnan");
md_smprod=mean(MD.prod,2,"omitnan");
lp_smprod=mean(LP.prod,2,"omitnan");
ld_smprod=mean(LD.prod,2,"omitnan");

sf_smnu=mean(SF.nu,2,"omitnan");
sp_smnu=mean(SP.nu,2,"omitnan");
sd_smnu=mean(SD.nu,2,"omitnan");
mf_smnu=mean(MF.nu,2,"omitnan");
mp_smnu=mean(MP.nu,2,"omitnan");
md_smnu=mean(MD.nu,2,"omitnan");
lp_smnu=mean(LP.nu,2,"omitnan");
ld_smnu=mean(LD.nu,2,"omitnan");

%% Each year
a = 1:12:nt; % start of each yr
b = 12:12:nt; % end of each yr
ympB = NaN*ones(length(sf_smprod),(nt/12));
ympSF = ympB;
ympSP = ympB;
ympSD = ympB;
ympMF = ympB;
ympMP = ympB;
ympMD = ympB;
ympLP = ympB;
ympLD = ympB;
ymnSF = ympB;
ymnSP = ympB;
ymnSD = ympB;
ymnMF = ympB;
ymnMP = ympB;
ymnMD = ympB;
ymnLP = ympB;
ymnLD = ympB;
for i = 1:(nt/12)
    %yr = (i+1987);
    %! Put vars of netcdf file
    ympSF(:,i) = mean(SF.prod(:,a(i):b(i)),2,"omitnan");
    ympSP(:,i) = mean(SP.prod(:,a(i):b(i)),2,"omitnan");
    ympSD(:,i) = mean(SD.prod(:,a(i):b(i)),2,"omitnan");
    ympMF(:,i) = mean(MF.prod(:,a(i):b(i)),2,"omitnan");
    ympMP(:,i) = mean(MP.prod(:,a(i):b(i)),2,"omitnan");
    ympMD(:,i) = mean(MD.prod(:,a(i):b(i)),2,"omitnan");
    ympLP(:,i) = mean(LP.prod(:,a(i):b(i)),2,"omitnan");
    ympLD(:,i) = mean(LD.prod(:,a(i):b(i)),2,"omitnan");

    ymnSF(:,i) = mean(SF.nu(:,a(i):b(i)),2,"omitnan");
    ymnSP(:,i) = mean(SP.nu(:,a(i):b(i)),2,"omitnan");
    ymnSD(:,i) = mean(SD.nu(:,a(i):b(i)),2,"omitnan");
    ymnMF(:,i) = mean(MF.nu(:,a(i):b(i)),2,"omitnan");
    ymnMP(:,i) = mean(MP.nu(:,a(i):b(i)),2,"omitnan");
    ymnMD(:,i) = mean(MD.nu(:,a(i):b(i)),2,"omitnan");
    ymnLP(:,i) = mean(LP.nu(:,a(i):b(i)),2,"omitnan");
    ymnLD(:,i) = mean(LD.nu(:,a(i):b(i)),2,"omitnan");
end

%%
save([fpath 'Means_Hist_' vers '_' harv '_prod_' cfile '.mat'],...
    'sf_smprod','sp_smprod','sd_smprod','mf_smprod','mp_smprod','md_smprod',...
    'lp_smprod','ld_smprod',...
    'sf_tmprod','sp_tmprod','sd_tmprod','mf_tmprod','mp_tmprod','md_tmprod',...
    'lp_tmprod','ld_tmprod',...
    'sf_tmnu','sp_tmnu','sd_tmnu','mf_tmnu','mp_tmnu','md_tmnu',...
    'lp_tmnu','ld_tmnu',...
    'sf_smnu','sp_smnu','sd_smnu','mf_smnu','mp_smnu','md_smnu',...
    'lp_smnu','ld_smnu',...
    'ympSF','ympSP','ympSD','ympMF','ympMP','ympMD','ympLP','ympLD',...
    'ymnSF','ymnSP','ymnSD','ymnMF','ymnMP','ymnMD','ymnLP','ymnLD');



