function netcdf_read_hist_fished_prod_ens(fname,simname)

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
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%Time
sp_tprod=mean(SP.prod,1);
sf_tprod=mean(SF.prod,1);
sd_tprod=mean(SD.prod,1);
mp_tprod=mean(MP.prod,1);
mf_tprod=mean(MF.prod,1);
md_tprod=mean(MD.prod,1);
lp_tprod=mean(LP.prod,1);
ld_tprod=mean(LD.prod,1);


%% 50 yrs (1951-2000)
y = 1860+(1/12):(1/12):2005;
yr50=find(y>=1951 & y<2001);
sp_prod50=nanmean(SP.prod(:,yr50),2);
sf_prod50=nanmean(SF.prod(:,yr50),2);
sd_prod50=nanmean(SD.prod(:,yr50),2);
mp_prod50=nanmean(MP.prod(:,yr50),2);
mf_prod50=nanmean(MF.prod(:,yr50),2);
md_prod50=nanmean(MD.prod(:,yr50),2);
lp_prod50=nanmean(LP.prod(:,yr50),2);
ld_prod50=nanmean(LD.prod(:,yr50),2);


% 1990-1995 (Climatology)
lyr=find(y>=1990 & y<1995);
%Means
sp_prod5=mean(SP.prod(:,lyr),2);
sf_prod5=mean(SF.prod(:,lyr),2);
sd_prod5=mean(SD.prod(:,lyr),2);
mp_prod5=mean(MP.prod(:,lyr),2);
mf_prod5=mean(MF.prod(:,lyr),2);
md_prod5=mean(MD.prod(:,lyr),2);
lp_prod5=mean(LP.prod(:,lyr),2);
ld_prod5=mean(LD.prod(:,lyr),2);


%% Every 5 years
st=1:60:length(time);
en=60:60:length(time);

for n=1:length(st)
    sp_prod(:,n)=nanmean(SP.prod(:,st(n):en(n)),2);
    sf_prod(:,n)=nanmean(SF.prod(:,st(n):en(n)),2);
    sd_prod(:,n)=nanmean(SD.prod(:,st(n):en(n)),2);
    mp_prod(:,n)=nanmean(MP.prod(:,st(n):en(n)),2);
    mf_prod(:,n)=nanmean(MF.prod(:,st(n):en(n)),2);
    md_prod(:,n)=nanmean(MD.prod(:,st(n):en(n)),2);
    lp_prod(:,n)=nanmean(LP.prod(:,st(n):en(n)),2);
    ld_prod(:,n)=nanmean(LD.prod(:,st(n):en(n)),2);
    
%     sp_rec(:,n)=nanmean(SP.rec(:,st(n):en(n)),2);
%     sf_rec(:,n)=nanmean(SF.rec(:,st(n):en(n)),2);
%     sd_rec(:,n)=nanmean(SD.rec(:,st(n):en(n)),2);
%     mp_rec(:,n)=nanmean(MP.rec(:,st(n):en(n)),2);
%     mf_rec(:,n)=nanmean(MF.rec(:,st(n):en(n)),2);
%     md_rec(:,n)=nanmean(MD.rec(:,st(n):en(n)),2);
%     lp_rec(:,n)=nanmean(LP.rec(:,st(n):en(n)),2);
%     ld_rec(:,n)=nanmean(LD.rec(:,st(n):en(n)),2);
    
%     mp_my(:,n)=nanmean(MP.yield(:,st(n):en(n)),2);
%     mf_my(:,n)=nanmean(MF.yield(:,st(n):en(n)),2);
%     md_my(:,n)=nanmean(MD.yield(:,st(n):en(n)),2);
%     lp_my(:,n)=nanmean(LP.yield(:,st(n):en(n)),2);
%     ld_my(:,n)=nanmean(LD.yield(:,st(n):en(n)),2);
end

%%
save([fname '_Means_prod_' simname '.mat'],'time',...
    'y','yr50','lyr',...
    'sf_tprod','sp_tprod','sd_tprod',...
    'mf_tprod','mp_tprod','md_tprod',...
    'lp_tprod','ld_tprod',...
    'sf_prod50','sp_prod50','sd_prod50',...
    'mf_prod50','mp_prod50','md_prod50',...
    'lp_prod50','ld_prod50',...
    'sf_prod5','sp_prod5','sd_prod5',...
    'mf_prod5','mp_prod5','md_prod5',...
    'lp_prod5','ld_prod5',...
    'sf_prod','sp_prod','sd_prod',...
    'mf_prod','mp_prod','md_prod',...
    'lp_prod','ld_prod');




end
