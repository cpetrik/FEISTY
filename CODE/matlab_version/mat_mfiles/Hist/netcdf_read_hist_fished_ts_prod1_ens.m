function netcdf_read_hist_fished_ts_prod1_ens(fname,simname)

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

% Benthic material (Time)
ncid = netcdf.open([fname '_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);
clear biomass

%% Take means and totals
y = 1860+(1/12):(1/12):2005;
yr50=find(y>=1951 & y<2001);
% 1990-1995 (Climatology)
lyr=find(y>=1990 & y<1995);

%% Every year
st=1:60:length(time);
en=60:60:length(time);

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
