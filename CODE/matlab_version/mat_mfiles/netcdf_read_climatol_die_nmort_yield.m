% POEM output at all locations

clear all
close all


cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'SWmlog_All_fish03';
vars = '_die_nmort_yield';

fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'Climatol_' harv vars '_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.mort = mort;
SP.die = die;
clear  yield die

% SF
ncid = netcdf.open([fpath 'Climatol_' harv vars '_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.mort = mort;
SF.die = die;
clear  yield die

% SD
ncid = netcdf.open([fpath 'Climatol_' harv vars '_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.mort = mort;
SD.die = die;
clear  yield die

% MP
ncid = netcdf.open([fpath 'Climatol_' harv vars '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.mort = mort;
MP.yield = yield;
MP.die = die;
clear  yield die

% MF
ncid = netcdf.open([fpath 'Climatol_' harv vars '_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.mort = mort;
MF.yield = yield;
MF.die = die;
clear  yield die mort

% MD
ncid = netcdf.open([fpath 'Climatol_' harv vars '_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.mort = mort;
MD.yield = yield;
MD.die = die;
clear  yield die

% LP
ncid = netcdf.open([fpath 'Climatol_' harv vars '_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.mort = mort;
LP.yield = yield;
LP.die = die;
clear  yield die mort

% LD
ncid = netcdf.open([fpath 'Climatol_' harv vars '_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.mort = mort;
LD.yield = yield;
LD.die = die;
clear  yield die mort


%% Take means
[ids,nt] = size(LD.die);

%Time
mp_tmyield=mean(MP.yield,1);
mf_tmyield=mean(MF.yield,1);
md_tmyield=mean(MD.yield,1);
lp_tmyield=mean(LP.yield,1);
ld_tmyield=mean(LD.yield,1);

sp_tmdie=mean(SP.die,1);
sf_tmdie=mean(SF.die,1);
sd_tmdie=mean(SD.die,1);
mp_tmdie=mean(MP.die,1);
mf_tmdie=mean(MF.die,1);
md_tmdie=mean(MD.die,1);
lp_tmdie=mean(LP.die,1);
ld_tmdie=mean(LD.die,1);

sp_tmmort=mean(SP.mort,1);
sf_tmmort=mean(SF.mort,1);
sd_tmmort=mean(SD.mort,1);
mp_tmmort=mean(MP.mort,1);
mf_tmmort=mean(MF.mort,1);
md_tmmort=mean(MD.mort,1);
lp_tmmort=mean(LP.mort,1);
ld_tmmort=mean(LD.mort,1);

%% Last year
time=1:nt;
lyr=time((end-12+1):end);

mp_myield=mean(MP.yield(:,lyr),2);
mf_myield=mean(MF.yield(:,lyr),2);
md_myield=mean(MD.yield(:,lyr),2);
lp_myield=mean(LP.yield(:,lyr),2);
ld_myield=mean(LD.yield(:,lyr),2);

sp_mdie=mean(SP.die(:,lyr),2);
sf_mdie=mean(SF.die(:,lyr),2);
sd_mdie=mean(SD.die(:,lyr),2);
mp_mdie=mean(MP.die(:,lyr),2);
mf_mdie=mean(MF.die(:,lyr),2);
md_mdie=mean(MD.die(:,lyr),2);
lp_mdie=mean(LP.die(:,lyr),2);
ld_mdie=mean(LD.die(:,lyr),2);

sp_mmort=mean(SP.mort(:,lyr),2);
sf_mmort=mean(SF.mort(:,lyr),2);
sd_mmort=mean(SD.mort(:,lyr),2);
mp_mmort=mean(MP.mort(:,lyr),2);
mf_mmort=mean(MF.mort(:,lyr),2);
md_mmort=mean(MD.mort(:,lyr),2);
lp_mmort=mean(LP.mort(:,lyr),2);
ld_mmort=mean(LD.mort(:,lyr),2);

%% Just last year
sp_die=SP.die(:,lyr);
sf_die=SF.die(:,lyr);
sd_die=SD.die(:,lyr);
mp_die=MP.die(:,lyr);
mf_die=MF.die(:,lyr);
md_die=MD.die(:,lyr);
lp_die=LP.die(:,lyr);
ld_die=LD.die(:,lyr);

sp_mort=SP.mort(:,lyr);
sf_mort=SF.mort(:,lyr);
sd_mort=SD.mort(:,lyr);
mp_mort=MP.mort(:,lyr);
mf_mort=MF.mort(:,lyr);
md_mort=MD.mort(:,lyr);
lp_mort=LP.mort(:,lyr);
ld_mort=LD.mort(:,lyr);

mp_yield=MP.yield(:,lyr);
mf_yield=MF.yield(:,lyr);
md_yield=MD.yield(:,lyr);
lp_yield=LP.yield(:,lyr);
ld_yield=LD.yield(:,lyr);

%%
% save([fpath 'Means_Climatol_' harv '_' cfile '.mat'],...
%     'mf_tmyield','mp_tmyield','md_tmyield','lp_tmyield','ld_tmyield',...
%     'sf_tmdie','sp_tmdie','sd_tmdie','mf_tmdie','mp_tmdie',...
%     'md_tmdie','lp_tmdie','ld_tmdie',...
%     'sf_tmmort','sp_tmmort','sd_tmmort','mf_tmmort','mp_tmmort',...
%     'md_tmmort','lp_tmmort','ld_tmmort',...
%     'mf_myield','mp_myield','md_myield','lp_myield','ld_myield',...
%     'sf_mdie','sp_mdie','sd_mdie','mf_mdie','mp_mdie',...
%     'md_mdie','lp_mdie','ld_mdie',...
%     'sf_mmort','sp_mmort','sd_mmort','mf_mmort','mp_mmort','md_mmort',...
%     'lp_mmort','ld_mmort',...
%     'time','lyr','sf_die','sp_die','sd_die',...
%     'mf_die','mp_die','md_die',...
%     'lp_die','ld_die','sf_mort','sp_mort','sd_mort',...
%     'mf_mort','mp_mort','md_mort',...
%     'lp_mort','ld_mort',...
%     'mf_yield','mp_yield','md_yield',...
%     'lp_yield','ld_yield','-append');

save([fpath 'Means' vars '_Climatol_' harv '_' cfile '.mat'],...
    'mf_tmyield','mp_tmyield','md_tmyield','lp_tmyield','ld_tmyield',...
    'sf_tmdie','sp_tmdie','sd_tmdie','mf_tmdie','mp_tmdie',...
    'md_tmdie','lp_tmdie','ld_tmdie',...
    'sf_tmmort','sp_tmmort','sd_tmmort','mf_tmmort','mp_tmmort',...
    'md_tmmort','lp_tmmort','ld_tmmort',...
    'mf_myield','mp_myield','md_myield','lp_myield','ld_myield',...
    'sf_mdie','sp_mdie','sd_mdie','mf_mdie','mp_mdie',...
    'md_mdie','lp_mdie','ld_mdie',...
    'sf_mmort','sp_mmort','sd_mmort','mf_mmort','mp_mmort','md_mmort',...
    'lp_mmort','ld_mmort',...
    'time','lyr','sf_die','sp_die','sd_die',...
    'mf_die','mp_die','md_die',...
    'lp_die','ld_die','sf_mort','sp_mort','sd_mort',...
    'mf_mort','mp_mort','md_mort',...
    'lp_mort','ld_mort',...
    'mf_yield','mp_yield','md_yield',...
    'lp_yield','ld_yield');







