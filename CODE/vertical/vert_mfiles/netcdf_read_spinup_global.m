% Spinup for 150 yrs saved every 5 yrs
% Save last month for initializing runs

clear 
close all

%%
cfile = 'NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/petrik-lab/Feisty/NC/MOM6-1D/Global/offline_feisty/' cfile '/'];

exper = 'Global_spinup_COBALTv3_halfdeg_HPcap';

%% SP
ncid = netcdf.open([fpath exper '_All_fish03_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

biomass(biomass>1e19) = nan;

[nid,nt] = size(biomass);
SP.bio = biomass;

clear biomass 

%% SF
ncid = netcdf.open([fpath exper '_All_fish03_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);
biomass(biomass>1e19) = nan;

SF.bio = biomass;
Sml_f.bio = biomass(:,nt);

clear biomass

% SD
ncid = netcdf.open([fpath exper '_All_fish03_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);
biomass(biomass>1e19) = nan;

SD.bio = biomass;

clear biomass 

%% MP
ncid = netcdf.open([fpath exper '_All_fish03_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);
biomass(biomass>1e19) = nan;

MP.bio = biomass;

clear biomass 

% MF
ncid = netcdf.open([fpath exper '_All_fish03_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);
biomass(biomass>1e19) = nan;

MF.bio = biomass;

clear biomass

% MD
ncid = netcdf.open([fpath exper '_All_fish03_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);
biomass(biomass>1e19) = nan;

MD.bio = biomass;

clear biomass 

% LP
ncid = netcdf.open([fpath exper '_All_fish03_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);
biomass(biomass>1e19) = nan;

LP.bio = biomass;

clear biomass 

% LD
ncid = netcdf.open([fpath exper '_All_fish03_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);
biomass(biomass>1e19) = nan;

LD.bio = biomass;

clear biomass 

%% Benthic material
ncid = netcdf.open([fpath exper '_All_fish03_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);
biomass(biomass>1e19) = nan;

Bent.bio = biomass;

clear biomass 

%% Take means
time = 1:nt;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];


%Time
sf_tmean=squeeze(mean(SF.bio,1,'omitnan'));
sp_tmean=squeeze(mean(SP.bio,1,'omitnan'));
sd_tmean=squeeze(mean(SD.bio,1,'omitnan'));
mf_tmean=squeeze(mean(MF.bio,1,'omitnan'));
mp_tmean=squeeze(mean(MP.bio,1,'omitnan'));
md_tmean=squeeze(mean(MD.bio,1,'omitnan'));
lp_tmean=squeeze(mean(LP.bio,1,'omitnan'));
ld_tmean=squeeze(mean(LD.bio,1,'omitnan'));
b_tmean=squeeze(mean(Bent.bio,1,'omitnan'));


%% Last year
% lyr=time((end-12+1):end);
% sp_mean=mean(SP.bio(:,lyr),2);
% sf_mean=mean(SF.bio(:,lyr),2);
% sd_mean=mean(SD.bio(:,lyr),2);
% mp_mean=mean(MP.bio(:,lyr),2);
% mf_mean=mean(MF.bio(:,lyr),2);
% md_mean=mean(MD.bio(:,lyr),2);
% lp_mean=mean(LP.bio(:,lyr),2);
% ld_mean=mean(LD.bio(:,lyr),2);
% b_mean=mean(Bent.bio(:,lyr),2);
% 
% %% Save last month for initializing hindcast runs
% 
% Sml_f.bio = mean(SF.bio(:,end),2,'omitnan');
% Sml_p.bio = mean(SP.bio(:,end),2,'omitnan');
% Sml_d.bio = mean(SD.bio(:,end),2,'omitnan');
% Med_f.bio = mean(MF.bio(:,end),2,'omitnan');
% Med_p.bio = mean(MP.bio(:,end),2,'omitnan');
% Med_d.bio = mean(MD.bio(:,end),2,'omitnan');
% Lrg_p.bio = mean(LP.bio(:,end),2,'omitnan');
% Lrg_d.bio = mean(LD.bio(:,end),2,'omitnan');
% BENT.bio  = mean(Bent.bio(:,end),2,'omitnan');
% 
% save([fpath 'Last_mo_' exper '_All_fish03_' cfile '.mat'],'Sml_f','Sml_p','Sml_d',...
%     'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')
% 
% 
% %% Save means
% save([fpath 'Means_' exper '_All_fish03_' cfile '.mat'],'time','lyr',...
%     'sf_mean','sp_mean','sd_mean',...
%     'mf_mean','mp_mean','md_mean',...
%     'b_mean','lp_mean','ld_mean',...
%     'sf_tmean','sp_tmean','sd_tmean',...
%     'mf_tmean','mp_tmean','md_tmean',...
%     'b_tmean','lp_tmean','ld_tmean');
% 
% 
%%
time = 5:5:150;

figure
plot(time,sf_tmean)
title('SF')

figure
plot(time,sp_tmean)
title('SP')

figure
plot(time,sd_tmean)
title('SD')

figure
plot(time,mf_tmean)
title('MF')

figure
plot(time,mp_tmean)
title('MP')

figure
plot(time,md_tmean)
title('MD')

figure
plot(time,lp_tmean)
title('LP')

figure
plot(time,ld_tmean)
title('LD')

figure
plot(time,b_tmean)
title('B')
