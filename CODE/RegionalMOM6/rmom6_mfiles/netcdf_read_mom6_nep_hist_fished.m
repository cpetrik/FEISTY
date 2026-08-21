% MOM6-NEP10 run

clear 
close all

%%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/project/Feisty/NC/Matlab_new_size/' cfile '/NEP10/'];

exper = 'NEP10_Hindcast1993_no_move_All_fish03';
%exper = 'NEP10_Hindcast1993_no_move_obsfish';

%% MP
ncid = netcdf.open([fpath exper '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.yield = yield;

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

MF.yield = yield;

clear yield

% MD
ncid = netcdf.open([fpath exper '_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.yield = yield;

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

LP.yield = yield;

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

LD.yield = yield;

clear biomass yield

%% Take means
nt = length(time);
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];


%Time
mp_yield_tmean=mean(MP.bio,1);
mf_yield_tmean=mean(MF.bio,1);
md_yield_tmean=mean(MD.bio,1);
lp_yield_tmean=mean(LP.bio,1);
ld_yield_tmean=mean(LD.bio,1);

%% Each year
st=1:12:length(time);
en=12:12:length(time);

for n=1:length(st)
    mp_yield_asum(:,n)=sum(MP.yield(:,st(n):en(n)),2,'omitnan');
    mf_yield_asum(:,n)=sum(MF.yield(:,st(n):en(n)),2,'omitnan');
    md_yield_asum(:,n)=sum(MD.yield(:,st(n):en(n)),2,'omitnan');
    lp_yield_asum(:,n)=sum(LP.yield(:,st(n):en(n)),2,'omitnan');
    ld_yield_asum(:,n)=sum(LD.yield(:,st(n):en(n)),2,'omitnan');
end

%% Whole time period mean
mf_yield_smean=mean(MF.yield,2,'omitnan');
mp_yield_smean=mean(MP.yield,2,'omitnan');
md_yield_smean=mean(MD.yield,2,'omitnan');
lp_yield_smean=mean(LP.yield,2,'omitnan');
ld_yield_smean=mean(LD.yield,2,'omitnan');

%% Whole time period each month
% mf_yield=MF.yield;
% mp_yield=MP.yield;
% md_yield=MD.yield;
% lp_yield=LP.yield;
% ld_yield=LD.yield;

%% Save means
save([fpath 'Means_yield_' exper '_' cfile '.mat'],'time',...
    'mf_yield_tmean','mp_yield_tmean','md_yield_tmean',...
    'lp_yield_tmean','ld_yield_tmean',...
    'mf_yield_asum','mp_yield_asum','md_yield_asum',...
    'lp_yield_asum','ld_yield_asum',...
    'mf_yield_smean','mp_yield_smean','md_yield_smean',...
    'lp_yield_smean','ld_yield_smean'); %,,...
    % 'mf_yield','mp_yield','md_yield',...
    % 'lp_yield','ld_yield');
    

