% POEM output at all locations

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];

%nt=12*95;

%% SP
ncid = netcdf.open([fpath 'Forecast_',harv,'_sml_p.nc'],'NC_NOWRITE');
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
% SP.rec = rec;
% SP.clev = clev;
% SP.con = con;
% SP.die = die;
% SP.gamma = gamma;
% SP.nu = nu;

clear biomass clev con die gamma nu rec time prod

% SF
ncid = netcdf.open([fpath 'Forecast_',harv,'_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
SF.prod = prod(:,1:nt);
% SF.rec = rec(:,1:nt);
% SF.clev = clev(:,1:nt);
% SF.con = con(:,1:nt);
% SF.die = die(:,1:nt);
% SF.gamma = gamma(:,1:nt);
% SF.nu = nu(:,1:nt);

clear biomass clev con die gamma nu rec time prod

% SD
ncid = netcdf.open([fpath 'Forecast_',harv,'_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
SD.prod = prod;
% SD.rec = rec;
% SD.clev = clev;
% SD.con = con;
% SD.die = die;
% SD.gamma = gamma;
% SD.nu = nu;

clear biomass clev con die gamma nu rec time prod

% MP
ncid = netcdf.open([fpath 'Forecast_',harv,'_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
MP.prod = prod;
MP.yield = yield;
% MP.rec = rec;
% MP.clev = clev;
% MP.con = con;
% MP.die = die;
% MP.gamma = gamma;
% MP.nu = nu;

clear biomass clev con die gamma nu rec time prod yield

% MF
ncid = netcdf.open([fpath 'Forecast_',harv,'_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
MF.prod = prod;
MF.yield = yield;
% MF.rec = rec;
% MF.rep = rep;
% MF.clev = clev;
% MF.con = con;
% MF.die = die;
% MF.gamma = gamma;
% MF.nu = nu;

clear biomass clev con die gamma nu rec rep time prod yield

% MD
ncid = netcdf.open([fpath 'Forecast_',harv,'_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
MD.prod = prod;
MD.yield = yield;
% MD.rec = rec;
% MD.clev = clev;
% MD.con = con;
% MD.die = die;
% MD.gamma = gamma;
% MD.nu = nu;

clear biomass clev con die gamma nu rec time prod yield

% LP
ncid = netcdf.open([fpath 'Forecast_',harv,'_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
LP.prod = prod;
LP.yield = yield;
% LP.rec = rec;
% LP.rep = rep;
% LP.clev = clev;
% LP.con = con;
% LP.die = die;
% LP.gamma = gamma;
% LP.nu = nu;

clear biomass clev con die gamma nu rec rep time prod yield

% LD
ncid = netcdf.open([fpath 'Forecast_',harv,'_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
LD.prod = prod;
LD.yield = yield;
% LD.rec = rec;
% LD.rep = rep;
% LD.clev = clev;
% LD.con = con;
% LD.die = die;
% LD.gamma = gamma;
% LD.nu = nu;

clear biomass clev con die gamma nu rec rep time prod yield

% Benthic material
ncid = netcdf.open([fpath 'Forecast_',harv,'_bent.nc'],'NC_NOWRITE');
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

% sp_tprod=mean(SP.prod,1);
% sf_tprod=mean(SF.prod,1);
% sd_tprod=mean(SD.prod,1);
% mp_tprod=mean(MP.prod,1);
% mf_tprod=mean(MF.prod,1);
% md_tprod=mean(MD.prod,1);
% lp_tprod=mean(LP.prod,1);
% ld_tprod=mean(LD.prod,1);
% 
% sp_trec=mean(SP.rec,1);
% sf_trec=mean(SF.rec,1);
% sd_trec=mean(SD.rec,1);
% mp_trec=mean(MP.rec,1);
% mf_trec=mean(MF.rec,1);
% md_trec=mean(MD.rec,1);
% lp_trec=mean(LP.rec,1);
% ld_trec=mean(LD.rec,1);

mf_tmy=mean(MF.yield,1);
mp_tmy=mean(MP.yield,1);
md_tmy=mean(MD.yield,1);
lp_tmy=mean(LP.yield,1);
ld_tmy=mean(LD.yield,1);

%% Every 5 years
st=1:60:length(time);
en=60:60:length(time);

for n=1:length(st)
    sp_mean(:,n)=nanmean(SP.bio(:,st(n):en(n)),2);
    sf_mean(:,n)=nanmean(SF.bio(:,st(n):en(n)),2);
    sd_mean(:,n)=nanmean(SD.bio(:,st(n):en(n)),2);
    mp_mean(:,n)=nanmean(MP.bio(:,st(n):en(n)),2);
    mf_mean(:,n)=nanmean(MF.bio(:,st(n):en(n)),2);
    md_mean(:,n)=nanmean(MD.bio(:,st(n):en(n)),2);
    lp_mean(:,n)=nanmean(LP.bio(:,st(n):en(n)),2);
    ld_mean(:,n)=nanmean(LD.bio(:,st(n):en(n)),2);
    b_mean(:,n)=nanmean(Bent.bio(:,st(n):en(n)),2);
    
%     sp_prod(:,n)=nanmean(SP.prod(:,st(n):en(n)),2);
%     sf_prod(:,n)=nanmean(SF.prod(:,st(n):en(n)),2);
%     sd_prod(:,n)=nanmean(SD.prod(:,st(n):en(n)),2);
%     mp_prod(:,n)=nanmean(MP.prod(:,st(n):en(n)),2);
%     mf_prod(:,n)=nanmean(MF.prod(:,st(n):en(n)),2);
%     md_prod(:,n)=nanmean(MD.prod(:,st(n):en(n)),2);
%     lp_prod(:,n)=nanmean(LP.prod(:,st(n):en(n)),2);
%     ld_prod(:,n)=nanmean(LD.prod(:,st(n):en(n)),2);
%     
%     sp_rec(:,n)=nanmean(SP.rec(:,st(n):en(n)),2);
%     sf_rec(:,n)=nanmean(SF.rec(:,st(n):en(n)),2);
%     sd_rec(:,n)=nanmean(SD.rec(:,st(n):en(n)),2);
%     mp_rec(:,n)=nanmean(MP.rec(:,st(n):en(n)),2);
%     mf_rec(:,n)=nanmean(MF.rec(:,st(n):en(n)),2);
%     md_rec(:,n)=nanmean(MD.rec(:,st(n):en(n)),2);
%     lp_rec(:,n)=nanmean(LP.rec(:,st(n):en(n)),2);
%     ld_rec(:,n)=nanmean(LD.rec(:,st(n):en(n)),2);
%     
    mp_my(:,n)=nanmean(MP.yield(:,st(n):en(n)),2);
    mf_my(:,n)=nanmean(MF.yield(:,st(n):en(n)),2);
    md_my(:,n)=nanmean(MD.yield(:,st(n):en(n)),2);
    lp_my(:,n)=nanmean(LP.yield(:,st(n):en(n)),2);
    ld_my(:,n)=nanmean(LD.yield(:,st(n):en(n)),2);
end


%% Last 50 years
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

%Means & medians
lyr=yr50;
all_mean1=mean(SP.bio(:,lyr),2)+mean(SF.bio(:,lyr),2)+mean(SD.bio(:,lyr),2)+...
    mean(MP.bio(:,lyr),2)+mean(MF.bio(:,lyr),2)+mean(MD.bio(:,lyr),2)+...
    mean(LP.bio(:,lyr),2)+mean(LD.bio(:,lyr),2);
all_mean2=mean((SP.bio(:,lyr)+SF.bio(:,lyr)+SD.bio(:,lyr)+...
    MP.bio(:,lyr)+MF.bio(:,lyr)+MD.bio(:,lyr)+...
    LP.bio(:,lyr)+LD.bio(:,lyr)),2);
all_median1=median(SP.bio(:,lyr),2)+median(SF.bio(:,lyr),2)+median(SD.bio(:,lyr),2)+...
    median(MP.bio(:,lyr),2)+median(MF.bio(:,lyr),2)+median(MD.bio(:,lyr),2)+...
    median(LP.bio(:,lyr),2)+median(LD.bio(:,lyr),2);
all_median2=median((SP.bio(:,lyr)+SF.bio(:,lyr)+SD.bio(:,lyr)+...
    MP.bio(:,lyr)+MF.bio(:,lyr)+MD.bio(:,lyr)+...
    LP.bio(:,lyr)+LD.bio(:,lyr)),2);

% sp_prod50=nanmean(SP.prod(:,yr50),2);
% sf_prod50=nanmean(SF.prod(:,yr50),2);
% sd_prod50=nanmean(SD.prod(:,yr50),2);
% mp_prod50=nanmean(MP.prod(:,yr50),2);
% mf_prod50=nanmean(MF.prod(:,yr50),2);
% md_prod50=nanmean(MD.prod(:,yr50),2);
% lp_prod50=nanmean(LP.prod(:,yr50),2);
% ld_prod50=nanmean(LD.prod(:,yr50),2);
% 
% sp_rec50=nanmean(SP.rec(:,yr50),2);
% sf_rec50=nanmean(SF.rec(:,yr50),2);
% sd_rec50=nanmean(SD.rec(:,yr50),2);
% mp_rec50=nanmean(MP.rec(:,yr50),2);
% mf_rec50=nanmean(MF.rec(:,yr50),2);
% md_rec50=nanmean(MD.rec(:,yr50),2);
% lp_rec50=nanmean(LP.rec(:,yr50),2);
% ld_rec50=nanmean(LD.rec(:,yr50),2);

mf_my50=mean(MF.yield(:,yr50),2);
mp_my50=mean(MP.yield(:,yr50),2);
md_my50=mean(MD.yield(:,yr50),2);
lp_my50=mean(LP.yield(:,yr50),2);
ld_my50=mean(LD.yield(:,yr50),2);

%% Save means
% save([fpath 'Means_fore_',harv,'_' cfile '.mat'],...
%     'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean','lp_mean','ld_mean',...
%     'sf_prod','sp_prod','sd_prod','mf_prod','mp_prod','md_prod','lp_prod','ld_prod',...
%     'sf_rec','sp_rec','sd_rec','mf_rec','mp_rec','md_rec','lp_rec','ld_rec',...
%     'sf_mean50','sp_mean50','sd_mean50','mf_mean50','mp_mean50','md_mean50','b_mean50','lp_mean50','ld_mean50',...
%     'sf_prod50','sp_prod50','sd_prod50','mf_prod50','mp_prod50','md_prod50','lp_prod50','ld_prod50',...
%     'sf_rec50','sp_rec50','sd_rec50','mf_rec50','mp_rec50','md_rec50','lp_rec50','ld_rec50',...
%     'sf_tmean','sp_tmean','sd_tmean','mf_tmean','mp_tmean','md_tmean','b_tmean','lp_tmean','ld_tmean',...
%     'sf_tprod','sp_tprod','sd_tprod','mf_tprod','mp_tprod','md_tprod','lp_tprod','ld_tprod',...
%     'sf_trec','sp_trec','sd_trec','mf_trec','mp_trec','md_trec','lp_trec','ld_trec','time','yr50',...
%     'mf_tmy','mp_tmy','md_tmy','lp_tmy','ld_tmy',...
%     'mf_my','mp_my','md_my','lp_my','ld_my',...
%     'mf_my50','mp_my50','md_my50','lp_my50','ld_my50');

save([fpath 'Means_fore_',harv,'_' cfile '.mat']);



