function netcdf_ts_areaint_hist_fished_bio_prod_ens(fname,simname,area)

% FEISTY output at all locations

%% SP
prod=[];
ncid = netcdf.open([fname '_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(biomass);
SP.bio = biomass;
SP.prod = prod;
clear biomass prod

% SF
prod=[];
ncid = netcdf.open([fname '_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
SF.prod = prod;
clear biomass prod

% SD
prod=[];
ncid = netcdf.open([fname '_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
SD.prod = prod;
clear biomass prod

% MP
prod=[];
ncid = netcdf.open([fname '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
MP.yield = yield;
MP.prod = prod;
clear biomass yield prod

% MF
prod=[];
ncid = netcdf.open([fname '_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
MF.yield = yield;
MF.prod = prod;
clear biomass yield prod

% MD
prod=[];
ncid = netcdf.open([fname '_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
MD.yield = yield;
MD.prod = prod;
clear biomass yield prod

% LP
prod=[];
ncid = netcdf.open([fname '_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
LP.yield = yield;
LP.prod = prod;
clear biomass yield prod

% LD
prod=[];
ncid = netcdf.open([fname '_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
LD.yield = yield;
LD.prod = prod;
clear biomass yield prod

% Benthic material
ncid = netcdf.open([fname '_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
clear biomass

%% Take means and totals
%Time
sp_tamean=mean(SP.bio.*area,1);
sf_tamean=mean(SF.bio.*area,1);
sd_tamean=mean(SD.bio.*area,1);
mp_tamean=mean(MP.bio.*area,1);
mf_tamean=mean(MF.bio.*area,1);
md_tamean=mean(MD.bio.*area,1);
lp_tamean=mean(LP.bio.*area,1);
ld_tamean=mean(LD.bio.*area,1);
b_tamean=mean(Bent.bio.*area,1);

F = SF.bio + MF.bio;
P = SP.bio + MP.bio + LP.bio;
D = SD.bio + MD.bio + LD.bio;

tPD = nansum(P.*area) ./ nansum((P.*area)+(D.*area));
tFD = nansum(F.*area) ./ nansum((F.*area)+(D.*area));
tPelD = nansum((P.*area)+(F.*area)) ./ ...
    nansum((P.*area)+(F.*area)+(D.*area));
tDPel = nansum(D.*area) ./ nansum((P.*area)+(F.*area)+(D.*area));
tDP = nansum(D.*area) ./ nansum((D.*area)+(P.*area));
tDF = nansum(D.*area) ./ nansum((D.*area)+(F.*area));
tFP = nansum(F.*area) ./ nansum((F.*area)+(P.*area));
tPF = nansum(P.*area) ./ nansum((P.*area)+(F.*area));

%%
save([fname '_Means_' simname '.mat'],...
    'sf_tamean','sp_tamean','sd_tamean',...
    'mf_tamean','mp_tamean','md_tamean',...
    'lp_tamean','ld_tamean','b_tamean',...
    'tPD','tFD','tPelD','tDPel','tDP','tDF','tFP','tPF','-append');


end

