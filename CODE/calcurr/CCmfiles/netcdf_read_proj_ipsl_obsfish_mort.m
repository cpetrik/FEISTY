% FEISTY output at all locations

clear
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';
vers = 'IPSL';
harv = 'All_fishobs';
%Project_IPSL_All_fishobs_bent.nc

fpath=['/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/',vers,'/'];

%% SP
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_mort_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[nid,nt] = size(die);

SP.die = die;
SP.pred = pred;
clear die pred

%% SF
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_mort_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.die = die(:,1:nt);
SF.pred = pred;
clear die pred

% SD
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_mort_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.die = die;
SD.pred = pred;
clear die pred

% MP
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_mort_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.die = die;
MP.pred = pred;
clear die pred

% MF
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_mort_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.die = die;
MF.pred = pred;
clear die pred

% MD
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_mort_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.die = die;
MD.pred = pred;
clear die pred

% LP
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_mort_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.die = die;
LP.pred = pred;
clear die pred

% LD
ncid = netcdf.open([fpath 'Project_' vers '_' harv '_mort_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.die = die;
LD.pred = pred;
clear die pred

%% Take means
%Time
sp_tmdie=mean(SP.die,1,"omitnan");
sf_tmdie=mean(SF.die,1,"omitnan");
sd_tmdie=mean(SD.die,1,"omitnan");
mp_tmdie=mean(MP.die,1,"omitnan");
mf_tmdie=mean(MF.die,1,"omitnan");
md_tmdie=mean(MD.die,1,"omitnan");
lp_tmdie=mean(LP.die,1,"omitnan");
ld_tmdie=mean(LD.die,1,"omitnan");

sf_tmpred=mean(SF.pred,1,"omitnan");
sp_tmpred=mean(SP.pred,1,"omitnan");
sd_tmpred=mean(SD.pred,1,"omitnan");
mf_tmpred=mean(MF.pred,1,"omitnan");
mp_tmpred=mean(MP.pred,1,"omitnan");
md_tmpred=mean(MD.pred,1,"omitnan");
lp_tmpred=mean(LP.pred,1,"omitnan");
ld_tmpred=mean(LD.pred,1,"omitnan");

%% All years
%lyr=time((end-12+1):end);
%lyr=1:12;
sp_smdie=mean(SP.die,2,"omitnan");
sf_smdie=mean(SF.die,2,"omitnan");
sd_smdie=mean(SD.die,2,"omitnan");
mp_smdie=mean(MP.die,2,"omitnan");
mf_smdie=mean(MF.die,2,"omitnan");
md_smdie=mean(MD.die,2,"omitnan");
lp_smdie=mean(LP.die,2,"omitnan");
ld_smdie=mean(LD.die,2,"omitnan");

sf_smpred=mean(SF.pred,2,"omitnan");
sp_smpred=mean(SP.pred,2,"omitnan");
sd_smpred=mean(SD.pred,2,"omitnan");
mf_smpred=mean(MF.pred,2,"omitnan");
mp_smpred=mean(MP.pred,2,"omitnan");
md_smpred=mean(MD.pred,2,"omitnan");
lp_smpred=mean(LP.pred,2,"omitnan");
ld_smpred=mean(LD.pred,2,"omitnan");

%% Each year
a = 1:12:nt; % start of each yr
b = 12:12:nt; % end of each yr
ymdB = NaN*ones(length(sf_smdie),(nt/12));
ymdSF = ymdB;
ymdSP = ymdB;
ymdSD = ymdB;
ymdMF = ymdB;
ymdMP = ymdB;
ymdMD = ymdB;
ymdLP = ymdB;
ymdLD = ymdB;
ymmSF = ymdB;
ymmSP = ymdB;
ymmSD = ymdB;
ymmMF = ymdB;
ymmMP = ymdB;
ymmMD = ymdB;
ymmLP = ymdB;
ymmLD = ymdB;
for i = 1:(nt/12)
    %yr = (i+1987);
    %! Put vars of netcdf file
    ymdSF(:,i) = mean(SF.die(:,a(i):b(i)),2,"omitnan");
    ymdSP(:,i) = mean(SP.die(:,a(i):b(i)),2,"omitnan");
    ymdSD(:,i) = mean(SD.die(:,a(i):b(i)),2,"omitnan");
    ymdMF(:,i) = mean(MF.die(:,a(i):b(i)),2,"omitnan");
    ymdMP(:,i) = mean(MP.die(:,a(i):b(i)),2,"omitnan");
    ymdMD(:,i) = mean(MD.die(:,a(i):b(i)),2,"omitnan");
    ymdLP(:,i) = mean(LP.die(:,a(i):b(i)),2,"omitnan");
    ymdLD(:,i) = mean(LD.die(:,a(i):b(i)),2,"omitnan");

    ymmSF(:,i) = mean(SF.pred(:,a(i):b(i)),2,"omitnan");
    ymmSP(:,i) = mean(SP.pred(:,a(i):b(i)),2,"omitnan");
    ymmSD(:,i) = mean(SD.pred(:,a(i):b(i)),2,"omitnan");
    ymmMF(:,i) = mean(MF.pred(:,a(i):b(i)),2,"omitnan");
    ymmMP(:,i) = mean(MP.pred(:,a(i):b(i)),2,"omitnan");
    ymmMD(:,i) = mean(MD.pred(:,a(i):b(i)),2,"omitnan");
    ymmLP(:,i) = mean(LP.pred(:,a(i):b(i)),2,"omitnan");
    ymmLD(:,i) = mean(LD.pred(:,a(i):b(i)),2,"omitnan");
end

%%
die_name = 'total biomass eaten';
pred_name = 'predation rate';
save([fpath 'Means_Project_' vers '_' harv '_mort_' cfile '.mat'],...
    'sf_smdie','sp_smdie','sd_smdie','mf_smdie','mp_smdie','md_smdie',...
    'lp_smdie','ld_smdie',...
    'sf_tmdie','sp_tmdie','sd_tmdie','mf_tmdie','mp_tmdie','md_tmdie',...
    'lp_tmdie','ld_tmdie',...
    'sf_tmpred','sp_tmpred','sd_tmpred','mf_tmpred','mp_tmpred','md_tmpred',...
    'lp_tmpred','ld_tmpred',...
    'sf_smpred','sp_smpred','sd_smpred','mf_smpred','mp_smpred','md_smpred',...
    'lp_smpred','ld_smpred',...
    'ymdSF','ymdSP','ymdSD','ymdMF','ymdMP','ymdMD','ymdLP','ymdLD',...
    'ymmSF','ymmSP','ymmSD','ymmMF','ymmMP','ymmMD','ymmLP','ymmLD',...
    'die_name','pred_name');



