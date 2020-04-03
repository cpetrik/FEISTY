function netcdf_read_fore_fished_nu_ens(fname,simname)

% FEISTY output at all locations

%% SP
ncid = netcdf.open([fname '_nu_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(nu);
SP.nu = nu;
clear nu

%% SF
ncid = netcdf.open([fname '_nu_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.nu = nu(:,1:nt);
clear nu

%% SD
ncid = netcdf.open([fname '_nu_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.nu = nu;
clear nu

% MP
ncid = netcdf.open([fname '_nu_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.nu = nu;
clear nu

% MF
ncid = netcdf.open([fname '_nu_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.nu = nu;
clear nu

% MD
ncid = netcdf.open([fname '_nu_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.nu = nu;
clear nu

% LP
ncid = netcdf.open([fname '_nu_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.nu = nu;
clear nu

% LD
ncid = netcdf.open([fname '_nu_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.nu = nu;
clear nu

% Benthic material
ncid = netcdf.open([fname '_nu_bent.nc'],'NC_NOWRITE');
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
sp_tnu=mean(SP.nu,1);
sf_tnu=mean(SF.nu,1);
sd_tnu=mean(SD.nu,1);
mp_tnu=mean(MP.nu,1);
mf_tnu=mean(MF.nu,1);
md_tnu=mean(MD.nu,1);
lp_tnu=mean(LP.nu,1);
ld_tnu=mean(LD.nu,1);

%% Last 50 years
yr50=time((end-(50*12)+1):end);
sp_nu50=mean(SP.nu(:,yr50),2);
sf_nu50=mean(SF.nu(:,yr50),2);
sd_nu50=mean(SD.nu(:,yr50),2);
mp_nu50=mean(MP.nu(:,yr50),2);
mf_nu50=mean(MF.nu(:,yr50),2);
md_nu50=mean(MD.nu(:,yr50),2);
lp_nu50=mean(LP.nu(:,yr50),2);
ld_nu50=mean(LD.nu(:,yr50),2);

%% Every 5 years
st=1:60:length(time);
en=60:60:length(time);

for n=1:length(st)
    sp_nu(:,n)=nanmean(SP.nu(:,st(n):en(n)),2);
    sf_nu(:,n)=nanmean(SF.nu(:,st(n):en(n)),2);
    sd_nu(:,n)=nanmean(SD.nu(:,st(n):en(n)),2);
    mp_nu(:,n)=nanmean(MP.nu(:,st(n):en(n)),2);
    mf_nu(:,n)=nanmean(MF.nu(:,st(n):en(n)),2);
    md_nu(:,n)=nanmean(MD.nu(:,st(n):en(n)),2);
    lp_nu(:,n)=nanmean(LP.nu(:,st(n):en(n)),2);
    ld_nu(:,n)=nanmean(LD.nu(:,st(n):en(n)),2);
    
end

%%
save([fname '_Means_nu_' simname '.mat'],'time','yr50',...
    'sf_tnu','sp_tnu','sd_tnu',...
    'mf_tnu','mp_tnu','md_tnu',...
    'lp_tnu','ld_tnu',...
    'sf_nu50','sp_nu50','sd_nu50',...
    'mf_nu50','mp_nu50','md_nu50',...
    'lp_nu50','ld_nu50',...
    'sf_nu','sp_nu','sd_nu',...
    'mf_nu','mp_nu','md_nu',...
    'lp_nu','ld_nu');

end

