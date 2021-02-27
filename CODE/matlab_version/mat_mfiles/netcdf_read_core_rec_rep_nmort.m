% POEM output at all locations

clear all
close all


cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
var = '_rec_rep_nmort';

fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'Core_' harv var '_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.mort = mort;
SP.rec = rec;
clear mort rec

% SF
ncid = netcdf.open([fpath 'Core_' harv var '_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.mort = mort;
SF.rec = rec;
clear mort rec

% SD
ncid = netcdf.open([fpath 'Core_' harv var '_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.mort = mort;
SD.rec = rec;
clear mort rec

% MP
ncid = netcdf.open([fpath 'Core_' harv var '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.mort = mort;
MP.rec = rec;
clear mort rec

% MF
ncid = netcdf.open([fpath 'Core_' harv var '_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.mort = mort;
MF.rep = rep;
MF.rec = rec;
clear mort rec rep

% MD
ncid = netcdf.open([fpath 'Core_' harv var '_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.mort = mort;
MD.rec = rec;
clear mort rec

% LP
ncid = netcdf.open([fpath 'Core_' harv var '_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.mort = mort;
LP.rep = rep;
LP.rec = rec;
clear mort rec rep

% LD
ncid = netcdf.open([fpath 'Core_' harv var '_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.mort = mort;
LD.rep = rep;
LD.rec = rec;
clear mort rec rep


%% Take means
[ids,nt] = size(LD.mort);

%Time
sp_tmmort=mean(SP.mort,1);
sf_tmmort=mean(SF.mort,1);
sd_tmmort=mean(SD.mort,1);
mp_tmmort=mean(MP.mort,1);
mf_tmmort=mean(MF.mort,1);
md_tmmort=mean(MD.mort,1);
lp_tmmort=mean(LP.mort,1);
ld_tmmort=mean(LD.mort,1);

sp_tmrec=mean(SP.rec,1);
sf_tmrec=mean(SF.rec,1);
sd_tmrec=mean(SD.rec,1);
mp_tmrec=mean(MP.rec,1);
mf_tmrec=mean(MF.rec,1);
md_tmrec=mean(MD.rec,1);
lp_tmrec=mean(LP.rec,1);
ld_tmrec=mean(LD.rec,1);

mf_tmrep=mean(MF.rep,1);
lp_tmrep=mean(LP.rep,1);
ld_tmrep=mean(LD.rep,1);

%% Spatial mean over all time
time=1:nt;
lyr=time((end-12+1):end);

sp_mmort=mean(SP.mort,2);
sf_mmort=mean(SF.mort,2);
sd_mmort=mean(SD.mort,2);
mp_mmort=mean(MP.mort,2);
mf_mmort=mean(MF.mort,2);
md_mmort=mean(MD.mort,2);
lp_mmort=mean(LP.mort,2);
ld_mmort=mean(LD.mort,2);

sp_mrec=mean(SP.rec,2);
sf_mrec=mean(SF.rec,2);
sd_mrec=mean(SD.rec,2);
mp_mrec=mean(MP.rec,2);
mf_mrec=mean(MF.rec,2);
md_mrec=mean(MD.rec,2);
lp_mrec=mean(LP.rec,2);
ld_mrec=mean(LD.rec,2);

mf_mrep=mean(MF.rep,2);
lp_mrep=mean(LP.rep,2);
ld_mrep=mean(LD.rep,2);

%% Every mo
sp_rec=SP.rec;
sf_rec=SF.rec;
sd_rec=SD.rec;
mp_rec=MP.rec;
mf_rec=MF.rec;
md_rec=MD.rec;
lp_rec=LP.rec;
ld_rec=LD.rec;

sp_mort=SP.mort;
sf_mort=SF.mort;
sd_mort=SD.mort;
mp_mort=MP.mort;
mf_mort=MF.mort;
md_mort=MD.mort;
lp_mort=LP.mort;
ld_mort=LD.mort;

mf_rep=MF.rep;
lp_rep=LP.rep;
ld_rep=LD.rep;

%%
% save([fpath 'Means_Core_' harv '_' cfile '.mat'],'sf_mmort',...
%     'sf_tmmort','sp_tmmort','sd_tmmort','mf_tmmort','mp_tmmort',...
%     'md_tmmort','lp_tmmort','ld_tmmort',...
%     'sf_tmrec','sp_tmrec','sd_tmrec','mf_tmrec','mp_tmrec',...
%     'md_tmrec','lp_tmrec','ld_tmrec',...
%     'mf_tmrep','lp_tmrep','ld_tmrep',...
%     'sf_mmort','sp_mmort','sd_mmort','mf_mmort','mp_mmort',...
%     'md_mmort','lp_mmort','ld_mmort',...
%     'sf_mrec','sp_mrec','sd_mrec','mf_mrec','mp_mrec',...
%     'md_mrec','lp_mrec','ld_mrec',...
%     'mf_mrep','lp_mrep','ld_mrep',...
%     'time','lyr','sf_rec','sp_rec','sd_rec',...
%     'mf_rec','mp_rec','md_rec',...
%     'lp_rec','ld_rec','b_rec',...
%     'mf_rep','lp_rep','ld_rep','-append');

save([fpath 'Means',var,'_Core_' harv '_' cfile '.mat'],'sf_tmmort',...
    'sf_tmmort','sp_tmmort','sd_tmmort','mf_tmmort','mp_tmmort',...
    'md_tmmort','lp_tmmort','ld_tmmort',...
    'sf_tmrec','sp_tmrec','sd_tmrec','mf_tmrec','mp_tmrec',...
    'md_tmrec','lp_tmrec','ld_tmrec',...
    'mf_tmrep','lp_tmrep','ld_tmrep',...
    'sf_mmort','sp_mmort','sd_mmort','mf_mmort','mp_mmort',...
    'md_mmort','lp_mmort','ld_mmort',...
    'sf_mrec','sp_mrec','sd_mrec','mf_mrec','mp_mrec',...
    'md_mrec','lp_mrec','ld_mrec',...
    'mf_mrep','lp_mrep','ld_mrep',...
    'time','lyr','sf_rec','sp_rec','sd_rec',...
    'mf_rec','mp_rec','md_rec',...
    'lp_rec','ld_rec',...
    'mf_rep','lp_rep','ld_rep','sf_mort','sp_mort','sd_mort',...
    'mf_mort','mp_mort','md_mort','lp_mort','ld_mort');







