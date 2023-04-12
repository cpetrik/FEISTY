% FEISTY output at all locations

clear 
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_CMIP6/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'Hist_empHP_prod_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%%
[ni,nt] = size(prod);

SP.prod = prod;
clear prod

%% SF
ncid = netcdf.open([fpath 'Hist_empHP_prod_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.prod = prod(:,1:nt);
clear prod 

% SD
ncid = netcdf.open([fpath 'Hist_empHP_prod_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.prod = prod;
clear prod 

% MP
ncid = netcdf.open([fpath 'Hist_empHP_prod_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.prod = prod;
clear prod

% MF
ncid = netcdf.open([fpath 'Hist_empHP_prod_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.prod = prod;
clear prod

% MD
ncid = netcdf.open([fpath 'Hist_empHP_prod_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.prod = prod;
clear prod

% LP
ncid = netcdf.open([fpath 'Hist_empHP_prod_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.prod = prod;
clear prod

% LD
ncid = netcdf.open([fpath 'Hist_empHP_prod_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.prod = prod;
clear prod

%% time
ncid = netcdf.open([fpath 'Hist_time.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%% Take means for my own visualization

%Time
sp_tprod=mean(SP.prod,1);
sf_tprod=mean(SF.prod,1);
sd_tprod=mean(SD.prod,1);
mp_tprod=mean(MP.prod,1);
mf_tprod=mean(MF.prod,1);
md_tprod=mean(MD.prod,1);
lp_tprod=mean(LP.prod,1);
ld_tprod=mean(LD.prod,1);

%% Space
t=time;
mo=t/12;
mo=mo+1950;
yrP=find(mo>2000 & mo<=2010); 

sp_mprod=mean(SP.prod(:,yrP),2);
sf_mprod=mean(SF.prod(:,yrP),2);
sd_mprod=mean(SD.prod(:,yrP),2);
mp_mprod=mean(MP.prod(:,yrP),2);
mf_mprod=mean(MF.prod(:,yrP),2);
md_mprod=mean(MD.prod(:,yrP),2);
lp_mprod=mean(LP.prod(:,yrP),2);
ld_mprod=mean(LD.prod(:,yrP),2);

save([fpath 'Means_Hist_empHP_2000-2010_' cfile '.mat'],...
    'sf_tprod','sp_tprod','sd_tprod',...
    'mf_tprod','mp_tprod','md_tprod',...
    'lp_tprod','ld_tprod',...
    'sf_mprod','sp_mprod','sd_mprod',...
    'mf_mprod','mp_mprod','md_mprod',...
    'lp_mprod','ld_mprod','-append')

figure
plot(mo,log10(lp_tprod),'b'); hold on;
plot(mo,log10(mf_tprod),'r'); hold on;
plot(mo,log10(ld_tprod),'k'); hold on;

