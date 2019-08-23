% POEM output at all locations

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];

%nt=12*95;

%% SP
ncid = netcdf.open([fpath 'Forecast_pristine_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(biomass);

SP.bio = biomass;
clear biomass prod

% SF
ncid = netcdf.open([fpath 'Forecast_pristine_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
clear biomass prod

% SD
ncid = netcdf.open([fpath 'Forecast_pristine_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
clear biomass prod

% MP
ncid = netcdf.open([fpath 'Forecast_pristine_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
clear biomass prod

% MF
ncid = netcdf.open([fpath 'Forecast_pristine_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
clear biomass prod

% MD
ncid = netcdf.open([fpath 'Forecast_pristine_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
clear biomass prod

% LP
ncid = netcdf.open([fpath 'Forecast_pristine_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
clear biomass prod

% LD
ncid = netcdf.open([fpath 'Forecast_pristine_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
clear biomass prod

% Benthic material
ncid = netcdf.open([fpath 'Forecast_pristine_bent.nc'],'NC_NOWRITE');
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
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%Time
sp_tmean=mean(SP.bio,1);
sf_tmean=mean(SF.bio,1);
sd_tmean=mean(SD.bio,1);
mp_tmean=mean(MP.bio,1);
mf_tmean=mean(MF.bio,1);
md_tmean=mean(MD.bio,1);
lp_tmean=mean(LP.bio,1);
ld_tmean=mean(LD.bio,1);
b_tmean=mean(Bent.bio,1);


%% Last 50 years
y = 2005+(1/12):(1/12):2100;
yr50=time((end-(50*12)+1):end);

sp_mean50=mean(SP.bio(:,yr50),2);
sf_mean50=mean(SF.bio(:,yr50),2);
sd_mean50=mean(SD.bio(:,yr50),2);
mp_mean50=mean(MP.bio(:,yr50),2);
mf_mean50=mean(MF.bio(:,yr50),2);
md_mean50=mean(MD.bio(:,yr50),2);
lp_mean50=mean(LP.bio(:,yr50),2);
ld_mean50=mean(LD.bio(:,yr50),2);
b_mean50=mean(Bent.bio(:,yr50),2);

%% Every 5 years
% st=1:60:length(time);
% en=60:60:length(time);
% 
% for n=1:length(st)
%     sp_mean5y(:,n)=nanmean(SP.bio(:,st(n):en(n)),2);
%     sf_mean5y(:,n)=nanmean(SF.bio(:,st(n):en(n)),2);
%     sd_mean5y(:,n)=nanmean(SD.bio(:,st(n):en(n)),2);
%     mp_mean5y(:,n)=nanmean(MP.bio(:,st(n):en(n)),2);
%     mf_mean5y(:,n)=nanmean(MF.bio(:,st(n):en(n)),2);
%     md_mean5y(:,n)=nanmean(MD.bio(:,st(n):en(n)),2);
%     lp_mean5y(:,n)=nanmean(LP.bio(:,st(n):en(n)),2);
%     ld_mean5y(:,n)=nanmean(LD.bio(:,st(n):en(n)),2);
%     b_mean5y(:,n)=nanmean(Bent.bio(:,st(n):en(n)),2);
%     
% end

%% Every year
st=1:12:length(time);
en=12:12:length(time);
for m=1:length(en)
    yr1 = st(m):en(m);
    sp_mean(:,m)=mean(SP.bio(:,yr1),2);
    sf_mean(:,m)=mean(SF.bio(:,yr1),2);
    sd_mean(:,m)=mean(SD.bio(:,yr1),2);
    mp_mean(:,m)=mean(MP.bio(:,yr1),2);
    mf_mean(:,m)=mean(MF.bio(:,yr1),2);
    md_mean(:,m)=mean(MD.bio(:,yr1),2);
    lp_mean(:,m)=mean(LP.bio(:,yr1),2);
    ld_mean(:,m)=mean(LD.bio(:,yr1),2);
    b_mean(:,m)=mean(Bent.bio(:,yr1),2);
end


%%
save([fpath 'Means_Forecast_pristine_bio_' cfile '.mat'],'time','y','yr50',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean',...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50','b_mean50',...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'b_mean','lp_mean','ld_mean');

% 'sf_mean5','sp_mean5','sd_mean5',...
%     'mf_mean5','mp_mean5','md_mean5',...
%     'lp_mean5','ld_mean5','b_mean5',...
    
