% POEM output at all locations


clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'Historic_' harv '_prod_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(prod);

SP.prod = prod;
Sml_p.prod = prod(:,nt);
clear prod

%% SF
ncid = netcdf.open([fpath 'Historic_' harv '_prod_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

SF.prod = prod(:,1:nt);
Sml_f.prod = prod(:,nt);
clear prod

%% SD
ncid = netcdf.open([fpath 'Historic_' harv '_prod_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

SD.prod = prod;
Sml_d.prod = prod(:,nt);
clear prod

% MP
ncid = netcdf.open([fpath 'Historic_' harv '_prod_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

MP.prod = prod;
MP.yield = yield;
Med_p.prod = prod(:,nt);
clear prod yield

% MF
ncid = netcdf.open([fpath 'Historic_' harv '_prod_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

MF.prod = prod;
MF.yield = yield;
Med_f.prod = prod(:,nt);
clear prod yield

% MD
ncid = netcdf.open([fpath 'Historic_' harv '_prod_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

MD.prod = prod;
MD.yield = yield;
Med_d.prod = prod(:,nt);
clear prod yield

% LP
ncid = netcdf.open([fpath 'Historic_' harv '_prod_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

LP.prod = prod;
LP.yield = yield;
Lrg_p.prod = prod(:,nt);
clear prod yield

% LD
ncid = netcdf.open([fpath 'Historic_' harv '_prod_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

LD.prod = prod;
LD.yield = yield;
Lrg_d.prod = prod(:,nt);
clear prod yield

% Benthic material
ncid = netcdf.open([fpath 'Historic_' harv '_prod_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
BENT.bio = biomass(:,nt);
clear biomass

%% Take means and totals
% Every year
st=1:12:length(time);
en=12:12:length(time);

for n=1:length(st)
    sp_prod1(:,n)=nanmean(SP.prod(:,st(n):en(n)),2);
    sf_prod1(:,n)=nanmean(SF.prod(:,st(n):en(n)),2);
    sd_prod1(:,n)=nanmean(SD.prod(:,st(n):en(n)),2);
    mp_prod1(:,n)=nanmean(MP.prod(:,st(n):en(n)),2);
    mf_prod1(:,n)=nanmean(MF.prod(:,st(n):en(n)),2);
    md_prod1(:,n)=nanmean(MD.prod(:,st(n):en(n)),2);
    lp_prod1(:,n)=nanmean(LP.prod(:,st(n):en(n)),2);
    ld_prod1(:,n)=nanmean(LD.prod(:,st(n):en(n)),2);
end

%%
save([fpath 'Means_Historic_' harv '_prod_' cfile '.mat'],...
'sf_prod1','sp_prod1','sd_prod1',...
'mf_prod1','mp_prod1','md_prod1',...
'lp_prod1','ld_prod1','-append');
