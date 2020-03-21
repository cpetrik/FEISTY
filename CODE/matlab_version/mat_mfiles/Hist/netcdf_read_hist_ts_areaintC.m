% POEM output at all locations

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
epath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat','AREA_OCN');
grid = csvread([cpath 'grid_csv.csv']);
area = AREA_OCN(grid(:,1));
area_mat = repmat(area,1,145*12);

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

% SD
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

SF.prod(SF.prod(:)<0)=0;
SP.prod(SP.prod(:)<0)=0;
SD.prod(SD.prod(:)<0)=0;
MF.prod(MF.prod(:)<0)=0;
MP.prod(MP.prod(:)<0)=0;
MD.prod(MD.prod(:)<0)=0;
LP.prod(LP.prod(:)<0)=0;
LD.prod(LD.prod(:)<0)=0;

%Time
sp_tamean=nansum(SP.prod.*area,1);
sf_tamean=nansum(SF.prod.*area,1);
sd_tamean=nansum(SD.prod.*area,1);
mp_tamean=nansum(MP.prod.*area,1);
mf_tamean=nansum(MF.prod.*area,1);
md_tamean=nansum(MD.prod.*area,1);
lp_tamean=nansum(LP.prod.*area,1);
ld_tamean=nansum(LD.prod.*area,1);
b_tamean=nansum(Bent.bio.*area,1);

F = SF.prod + MF.prod;
P = SP.prod + MP.prod + LP.prod;
D = SD.prod + MD.prod + LD.prod;
All = F+P+D;

tPD = nansum(P.*area) ./ nansum((P.*area)+(D.*area));
tFD = nansum(F.*area) ./ nansum((F.*area)+(D.*area));
tPelD = nansum((P.*area)+(F.*area)) ./ ...
    nansum((P.*area)+(F.*area)+(D.*area));
tDPel = nansum(D.*area) ./ nansum((P.*area)+(F.*area)+(D.*area));
tDP = nansum(D.*area) ./ nansum((D.*area)+(P.*area));
tDF = nansum(D.*area) ./ nansum((D.*area)+(F.*area));
tFP = nansum(F.*area) ./ nansum((F.*area)+(P.*area));
tPF = nansum(P.*area) ./ nansum((P.*area)+(F.*area));

tP = nansum(P.*area) ./ nansum(All.*area);
tF = nansum(F.*area) ./ nansum(All.*area);
tPel = nansum((P.*area)+(F.*area)) ./ nansum(All.*area);

%%
save([fpath 'Means_Historic_' harv '_prod_' cfile '.mat'],...
    'sf_tamean','sp_tamean','sd_tamean',...
    'mf_tamean','mp_tamean','md_tamean',...
    'lp_tamean','ld_tamean','b_tamean',...
    'tPD','tFD','tPelD','tDPel','tDP','tDF','tFP','tPF',...
    'tP','tF','tPel','-append');

save([epath 'Means_Historic_' harv '_prod_' cfile '.mat'],...
    'sf_tamean','sp_tamean','sd_tamean',...
    'mf_tamean','mp_tamean','md_tamean',...
    'lp_tamean','ld_tamean','b_tamean',...
    'tPD','tFD','tPelD','tDPel','tDP','tDF','tFP','tPF',...
    'tP','tF','tPel','-append');

