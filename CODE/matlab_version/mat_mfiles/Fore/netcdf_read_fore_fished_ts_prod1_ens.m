function netcdf_read_fore_fished_ts_prod1_ens(fname,simname)

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

%% Take means 
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
save([fname '_Means_prod_' simname '.mat'],...
    'sf_prod1','sp_prod1','sd_prod1',...
    'mf_prod1','mp_prod1','md_prod1',...
    'lp_prod1','ld_prod1','-append');

end

