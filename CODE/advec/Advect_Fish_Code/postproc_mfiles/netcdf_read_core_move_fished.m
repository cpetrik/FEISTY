% CORE-forced run

clear 
close all

%%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/CORE/'];

exper = 'CORE_Hindcast_move_nu_v21_dt6h';

%% SP
ncid = netcdf.open([fpath exper '_All_fish03_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[nid,nt] = size(biomass);
SP.bio = biomass;

clear biomass 

%% SF
ncid = netcdf.open([fpath exper '_All_fish03_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);

clear biomass

% SD
ncid = netcdf.open([fpath exper '_All_fish03_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;

clear biomass 

%% MP
ncid = netcdf.open([fpath exper '_All_fish03_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
% MP.yield = yield;

clear biomass yield

% MF
ncid = netcdf.open([fpath exper '_All_fish03_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
% MF.yield = yield;

clear biomassyield

% MD
ncid = netcdf.open([fpath exper '_All_fish03_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
% MD.yield = yield;

clear biomass yield

% LP
ncid = netcdf.open([fpath exper '_All_fish03_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
% LP.yield = yield;

clear biomass yield

% LD
ncid = netcdf.open([fpath exper '_All_fish03_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
% LD.yield = yield;

clear biomass yield

% Benthic material
ncid = netcdf.open([fpath exper '_All_fish03_bent.nc'],'NC_NOWRITE');
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
sp_tmean=mean(SP.bio,1,'omitnan');
sf_tmean=mean(SF.bio,1,'omitnan');
sd_tmean=mean(SD.bio,1,'omitnan');
mp_tmean=mean(MP.bio,1,'omitnan');
mf_tmean=mean(MF.bio,1,'omitnan');
md_tmean=mean(MD.bio,1,'omitnan');
lp_tmean=mean(LP.bio,1,'omitnan');
ld_tmean=mean(LD.bio,1,'omitnan');
b_tmean=mean(Bent.bio,1,'omitnan');

%% Each year
st=1:12:length(time);
en=12:12:length(time);

for n=1:length(st)
    sp_smean(:,n)=mean(SP.bio(:,st(n):en(n)),2,'omitnan');
    sf_smean(:,n)=mean(SF.bio(:,st(n):en(n)),2,'omitnan');
    sd_smean(:,n)=mean(SD.bio(:,st(n):en(n)),2,'omitnan');
    mp_smean(:,n)=mean(MP.bio(:,st(n):en(n)),2,'omitnan');
    mf_smean(:,n)=mean(MF.bio(:,st(n):en(n)),2,'omitnan');
    md_smean(:,n)=mean(MD.bio(:,st(n):en(n)),2,'omitnan');
    lp_smean(:,n)=mean(LP.bio(:,st(n):en(n)),2,'omitnan');
    ld_smean(:,n)=mean(LD.bio(:,st(n):en(n)),2,'omitnan');
    b_smean(:,n)=mean(Bent.bio(:,st(n):en(n)),2,'omitnan');
    
%     mp_my(:,n)=mean(MP.yield(:,st(n):en(n)),2,'omitnan');
%     mf_my(:,n)=mean(MF.yield(:,st(n):en(n)),2,'omitnan');
%     md_my(:,n)=mean(MD.yield(:,st(n):en(n)),2,'omitnan');
%     lp_my(:,n)=mean(LP.yield(:,st(n):en(n)),2,'omitnan');
%     ld_my(:,n)=mean(LD.yield(:,st(n):en(n)),2,'omitnan');
end

%% Whole time period mean
sp_mean20=mean(SP.bio,2,'omitnan');
sf_mean20=mean(SF.bio,2,'omitnan');
sd_mean20=mean(SD.bio,2,'omitnan');
mp_mean20=mean(MP.bio,2,'omitnan');
mf_mean20=mean(MF.bio,2,'omitnan');
md_mean20=mean(MD.bio,2,'omitnan');
lp_mean20=mean(LP.bio,2,'omitnan');
ld_mean20=mean(LD.bio,2,'omitnan');
b_mean20=mean(Bent.bio,2,'omitnan');

% mf_my20=mean(MF.yield,2);
% mp_my20=mean(MP.yield,2);
% md_my20=mean(MD.yield,2);
% lp_my20=mean(LP.yield,2);
% ld_my20=mean(LD.yield,2);

%% Whole time period each month
% sp_bio=SP.bio;
% sf_bio=SF.bio;
% sd_bio=SD.bio;
% mp_bio=MP.bio;
% mf_bio=MF.bio;
% md_bio=MD.bio;
% lp_bio=LP.bio;
% ld_bio=LD.bio;
% b_bio=Bent.bio;
% 
% mf_yield=MF.yield;
% mp_yield=MP.yield;
% md_yield=MD.yield;
% lp_yield=LP.yield;
% ld_yield=LD.yield;

%% Save means
save([fpath 'Means_' exper '_' cfile '.mat'],'time',...
    'sf_smean','sp_smean','sd_smean',...
    'mf_smean','mp_smean','md_smean',...
    'b_smean','lp_smean','ld_smean',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'b_tmean','lp_tmean','ld_tmean',...
    'sf_mean20','sp_mean20','sd_mean20',...
    'mf_mean20','mp_mean20','md_mean20',...
    'lp_mean20','ld_mean20','b_mean20'); %,...
    % 'mf_my20','mp_my20','md_my20',...
    % 'lp_my20','ld_my20',...
    % 'sf_bio','sp_bio','sd_bio',...
    % 'mf_bio','mp_bio','md_bio',...
    % 'b_bio','lp_bio','ld_bio',...
    % 'mf_yield','mp_yield','md_yield',...
    % 'lp_yield','ld_yield');
    

