% POEM output at all locations

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'Historic_pristine_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(biomass);

SP.bio = biomass;
Sml_p.bio = biomass(:,nt);
clear biomass prod

% SF
ncid = netcdf.open([fpath 'Historic_pristine_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
Sml_f.bio = biomass(:,nt);
clear biomass prod

% SD
ncid = netcdf.open([fpath 'Historic_pristine_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
Sml_d.bio = biomass(:,nt);
clear biomass prod

% MP
ncid = netcdf.open([fpath 'Historic_pristine_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
Med_p.bio = biomass(:,nt);
clear biomass

% MF
ncid = netcdf.open([fpath 'Historic_pristine_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
Med_f.bio = biomass(:,nt);
clear biomass

% MD
ncid = netcdf.open([fpath 'Historic_pristine_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
Med_d.bio = biomass(:,nt);
clear biomass

% LP
ncid = netcdf.open([fpath 'Historic_pristine_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
Lrg_p.bio = biomass(:,nt);
clear biomass

% LD
ncid = netcdf.open([fpath 'Historic_pristine_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
Lrg_d.bio = biomass(:,nt);
clear biomass

% Benthic material
ncid = netcdf.open([fpath 'Historic_pristine_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
BENT.bio = biomass(:,nt);
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


%% 50 yrs (1951-2000)
y = 1860+(1/12):(1/12):2005;
yr50=find(y>=1951 & y<2001);
sp_mean50=mean(SP.bio(:,yr50),2);
sf_mean50=mean(SF.bio(:,yr50),2);
sd_mean50=mean(SD.bio(:,yr50),2);
mp_mean50=mean(MP.bio(:,yr50),2);
mf_mean50=mean(MF.bio(:,yr50),2);
md_mean50=mean(MD.bio(:,yr50),2);
lp_mean50=mean(LP.bio(:,yr50),2);
ld_mean50=mean(LD.bio(:,yr50),2);
b_mean50=mean(Bent.bio(:,yr50),2);

% 1990-1995 (Climatology)
lyr=find(y>=1990 & y<1995);
%Means
sp_mean5=mean(SP.bio(:,lyr),2);
sf_mean5=mean(SF.bio(:,lyr),2);
sd_mean5=mean(SD.bio(:,lyr),2);
mp_mean5=mean(MP.bio(:,lyr),2);
mf_mean5=mean(MF.bio(:,lyr),2);
md_mean5=mean(MD.bio(:,lyr),2);
lp_mean5=mean(LP.bio(:,lyr),2);
ld_mean5=mean(LD.bio(:,lyr),2);
b_mean5=mean(Bent.bio(:,lyr),2);
all_mean1=mean(SP.bio(:,lyr),2)+mean(SF.bio(:,lyr),2)+mean(SD.bio(:,lyr),2)+...
    mean(MP.bio(:,lyr),2)+mean(MF.bio(:,lyr),2)+mean(MD.bio(:,lyr),2)+...
    mean(LP.bio(:,lyr),2)+mean(LD.bio(:,lyr),2);
all_mean2=mean((SP.bio(:,lyr)+SF.bio(:,lyr)+SD.bio(:,lyr)+...
    MP.bio(:,lyr)+MF.bio(:,lyr)+MD.bio(:,lyr)+...
    LP.bio(:,lyr)+LD.bio(:,lyr)),2);
%Medians
all_median1=median(SP.bio(:,lyr),2)+median(SF.bio(:,lyr),2)+median(SD.bio(:,lyr),2)+...
    median(MP.bio(:,lyr),2)+median(MF.bio(:,lyr),2)+median(MD.bio(:,lyr),2)+...
    median(LP.bio(:,lyr),2)+median(LD.bio(:,lyr),2);
all_median2=median((SP.bio(:,lyr)+SF.bio(:,lyr)+SD.bio(:,lyr)+...
    MP.bio(:,lyr)+MF.bio(:,lyr)+MD.bio(:,lyr)+...
    LP.bio(:,lyr)+LD.bio(:,lyr)),2);

% 1990 (Climatology)
yr1=find(y>=1990 & y<1991);
%Totals
sp_tot=sum(SP.bio(:,yr1).*MNTH,2);
sf_tot=sum(SF.bio(:,yr1).*MNTH,2);
sd_tot=sum(SD.bio(:,yr1).*MNTH,2);
mp_tot=sum(MP.bio(:,yr1).*MNTH,2);
mf_tot=sum(MF.bio(:,yr1).*MNTH,2);
md_tot=sum(MD.bio(:,yr1).*MNTH,2);
lp_tot=sum(LP.bio(:,yr1).*MNTH,2);
ld_tot=sum(LD.bio(:,yr1).*MNTH,2);
b_tot=sum(Bent.bio(:,yr1).*MNTH,2);

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
    
end


%%
save([fpath 'Means_Historic_pristine_' cfile '.mat'],'time','y','yr50','yr1','lyr',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean',...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50','b_mean50',...
    'sf_mean5','sp_mean5','sd_mean5',...
    'mf_mean5','mp_mean5','md_mean5',...
    'lp_mean5','ld_mean5','b_mean5',...
    'sf_tot','sp_tot','sd_tot',...
    'mf_tot','mp_tot','md_tot',...
    'lp_tot','ld_tot','b_tot',...
    'all_median1','all_median2','all_mean1','all_mean2',...
    'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
    'lp_mean','ld_mean');


% Save last year for initializing forecast runs
save([fpath 'Last_mo_hist_pristine_' cfile '.mat'],'Sml_f','Sml_p','Sml_d',... 
    'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')





