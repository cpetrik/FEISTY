function netcdf_ts_areaint_hist_fished_prod_ens(fname,simname,area)

% FEISTY output at all locations

%% SP
prod=[];
ncid = netcdf.open([fname '_prod_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(prod);

SP.prod = prod;
clear prod

%% SF
prod=[];
ncid = netcdf.open([fname '_prod_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

SF.prod = prod(:,1:nt);
clear prod

% SD
prod=[];
ncid = netcdf.open([fname '_prod_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

SD.prod = prod;
clear prod

%% MP
prod=[];
ncid = netcdf.open([fname '_prod_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

MP.prod = prod;
clear prod yield

% MF
prod=[];
ncid = netcdf.open([fname '_prod_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

MF.prod = prod;
clear prod yield

% MD
prod=[];
ncid = netcdf.open([fname '_prod_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

MD.prod = prod;
clear prod yield

% LP
prod=[];
ncid = netcdf.open([fname '_prod_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

LP.prod = prod;
clear prod yield

% LD
prod=[];
ncid = netcdf.open([fname '_prod_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

LD.prod = prod;
clear prod yield

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

% Get rid of negative values
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
save([fname '_Means_prod_' simname '.mat'],...
    'sf_tamean','sp_tamean','sd_tamean',...
    'mf_tamean','mp_tamean','md_tamean',...
    'lp_tamean','ld_tamean','b_tamean',...
    'tPD','tFD','tPelD','tDPel','tDP','tDF','tFP','tPF',...
    'tP','tF','tPel','-append');


end

