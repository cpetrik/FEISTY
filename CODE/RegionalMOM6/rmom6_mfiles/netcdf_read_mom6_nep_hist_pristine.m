% MOM6-NEP10 run

clear 
close all

%%
%cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_CC80_RE00100';
fpath=['/project/Feisty/NC/Matlab_new_size/' cfile '/NEP10/'];

%exper = 'NEP10_Hindcast1993_no_move_All_fish03';
exper = 'NEP10_Hindcast1993_no_move_pristine';

%% SP
ncid = netcdf.open([fpath exper '_sml_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath exper '_sml_f.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath exper '_sml_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath exper '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
%MP.yield = yield;

clear biomass yield

% MF
ncid = netcdf.open([fpath exper '_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
%MF.yield = yield;

clear biomassyield

% MD
ncid = netcdf.open([fpath exper '_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
%MD.yield = yield;

clear biomass yield

% LP
ncid = netcdf.open([fpath exper '_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
%LP.yield = yield;

clear biomass yield

% LD
ncid = netcdf.open([fpath exper '_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
%LD.yield = yield;

clear biomass yield

% Benthic material
ncid = netcdf.open([fpath exper '_bent.nc'],'NC_NOWRITE');
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
    
    % mp_my(:,n)=mean(MP.yield(:,st(n):en(n)),2,'omitnan');
    % mf_my(:,n)=mean(MF.yield(:,st(n):en(n)),2,'omitnan');
    % md_my(:,n)=mean(MD.yield(:,st(n):en(n)),2,'omitnan');
    % lp_my(:,n)=mean(LP.yield(:,st(n):en(n)),2,'omitnan');
    % ld_my(:,n)=mean(LD.yield(:,st(n):en(n)),2,'omitnan');
end

%% Whole time period mean
sp_mean=mean(SP.bio,2,'omitnan');
sf_mean=mean(SF.bio,2,'omitnan');
sd_mean=mean(SD.bio,2,'omitnan');
mp_mean=mean(MP.bio,2,'omitnan');
mf_mean=mean(MF.bio,2,'omitnan');
md_mean=mean(MD.bio,2,'omitnan');
lp_mean=mean(LP.bio,2,'omitnan');
ld_mean=mean(LD.bio,2,'omitnan');
b_mean=mean(Bent.bio,2,'omitnan');

% mf_my=mean(MF.yield,2,'omitnan');
% mp_my=mean(MP.yield,2,'omitnan');
% md_my=mean(MD.yield,2,'omitnan');
% lp_my=mean(LP.yield,2,'omitnan');
% ld_my=mean(LD.yield,2,'omitnan');

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
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean'); %,...
    % 'mf_my','mp_my','md_my',...
    % 'lp_my','ld_my',,...
    % 'sf_bio','sp_bio','sd_bio',...
    % 'mf_bio','mp_bio','md_bio',...
    % 'b_bio','lp_bio','ld_bio',...
    % 'mf_yield','mp_yield','md_yield',...
    % 'lp_yield','ld_yield');
    

