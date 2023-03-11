% FEISTY output at all locations

clear all
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];

mod = 'v3.2';

%% time
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_time.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

st=1:12:length(time);
en=12:12:length(time);

%% SP
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(biomass);

SP.bio = biomass;
clear biomass

%Time
sp_tmean=mean(SP.bio,1,'omitnan');
%Space
sp_smean=mean(SP.bio,2,'omitnan');
% Each year
for n=1:length(st)
    sp_mean(:,n)=nanmean(SP.bio(:,st(n):en(n)),2);
end

clear SP

%% SF 
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
clear biomass 

%Time
sf_tmean=mean(SF.bio,1,'omitnan');
%Space
sf_smean=mean(SF.bio,2,'omitnan');
% Each year
for n=1:length(st)
    sf_mean(:,n)=nanmean(SF.bio(:,st(n):en(n)),2);
end

clear SF

%% SD
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
clear biomass 

%Time
sd_tmean=mean(SD.bio,1,'omitnan');
%Space
sd_smean=mean(SD.bio,2,'omitnan');
% Each year
for n=1:length(st)
    sd_mean(:,n)=nanmean(SD.bio(:,st(n):en(n)),2);
end

clear SD

%% MP
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
clear biomass

%Time
mp_tmean=mean(MP.bio,1,'omitnan');
%Space
mp_smean=mean(MP.bio,2,'omitnan');
% Each year
for n=1:length(st)
    mp_mean(:,n)=nanmean(MP.bio(:,st(n):en(n)),2);
end

clear MP

%% MF
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
clear biomass

%Time
mf_tmean=mean(MF.bio,1,'omitnan');
%Space
mf_smean=mean(MF.bio,2,'omitnan');
% Each year
for n=1:length(st)
    mf_mean(:,n)=nanmean(MF.bio(:,st(n):en(n)),2);
end

%clear MF

%% MD
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
clear biomass

%Time
md_tmean=mean(MD.bio,1,'omitnan');
%Space
md_smean=mean(MD.bio,2,'omitnan');
% Each year
for n=1:length(st)
    md_mean(:,n)=nanmean(MD.bio(:,st(n):en(n)),2);
end

clear MD

%% LP
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
clear biomass

%Time
lp_tmean=mean(LP.bio,1,'omitnan');
lp_smean=mean(LP.bio,2,'omitnan');
for n=1:length(st)
    lp_mean(:,n)=nanmean(LP.bio(:,st(n):en(n)),2);
end

clear LP

%% LD
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
clear biomass

%Time
ld_tmean=mean(LD.bio,1,'omitnan');
ld_smean=mean(LD.bio,2,'omitnan');
for n=1:length(st)
    ld_mean(:,n)=nanmean(LD.bio(:,st(n):en(n)),2);
end

clear LD

%% Benthic material
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
clear biomass

%Time
b_tmean=mean(Bent.bio,1,'omitnan');
b_smean =mean(Bent.bio,2,'omitnan');
for n=1:length(st)
    b_mean(:,n)=nanmean(Bent.bio(:,st(n):en(n)),2);
end

clear Bent

%% MF catch
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_catch_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.yield = yield;
clear yield

%Time
mf_tmcatch=mean(MF.yield,1,'omitnan');
mf_mcatch=mean(MF.yield,2,'omitnan');
for n=1:length(st)
    mf_my(:,n)=nanmean(MF.yield(:,st(n):en(n)),2);
end

clear MF

% MP catch
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_catch_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.yield = yield;
clear yield

mp_tmcatch=mean(MP.yield,1,'omitnan');
mp_mcatch=mean(MP.yield,2,'omitnan');
for n=1:length(st)
    mp_my(:,n)=nanmean(MP.yield(:,st(n):en(n)),2);
end

clear MP

% MD catch
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_catch_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.yield = yield;
clear yield

%Time
md_tmcatch=mean(MD.yield,1,'omitnan');
md_mcatch=mean(MD.yield,2,'omitnan');
for n=1:length(st)
    md_my(:,n)=nanmean(MD.yield(:,st(n):en(n)),2);
end

clear MD

% LP catch
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_catch_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.yield = yield;
clear yield

%Time
lp_tmcatch=mean(LP.yield,1,'omitnan');
lp_mcatch=mean(LP.yield,2,'omitnan');
for n=1:length(st)
    lp_my(:,n)=nanmean(LP.yield(:,st(n):en(n)),2);
end

clear LP

% LD catch
ncid = netcdf.open([fpath 'PI_ctrlclim_All_fishobs_',mod,'_empHP_catch_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.yield = yield;
clear yield

%Time
ld_tmcatch=mean(LD.yield,1,'omitnan');
ld_mcatch=mean(LD.yield,2,'omitnan');
for n=1:length(st)
    ld_my(:,n)=nanmean(LD.yield(:,st(n):en(n)),2);
end

clear LD

%
save([fpath 'Means_PI_ctrlclim_All_fishobs_',mod,'_' cfile '.mat'],'time',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean',...
    'mf_tmcatch','mp_tmcatch','md_tmcatch',...
    'lp_tmcatch','ld_tmcatch',...
    'sf_smean','sp_smean','sd_smean',...
    'mf_smean','mp_smean','md_smean',...
    'lp_smean','ld_smean','b_smean',...
    'mf_mcatch','mp_mcatch','md_mcatch',...
    'lp_mcatch','ld_mcatch',...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean',...
    'mf_my','mp_my','md_my',...
    'lp_my','ld_my')

%'mf_tmean','mf_smean','mf_mean',


%% Save last year for initializing forecast runs
Sml_f.bio = sf_mean(:,end);
Sml_p.bio = sp_mean(:,end);
Sml_d.bio = sd_mean(:,end);
Med_f.bio = mf_mean(:,end);
Med_p.bio = mp_mean(:,end);
Med_d.bio = md_mean(:,end);
Lrg_p.bio = lp_mean(:,end);
Lrg_d.bio = ld_mean(:,end);
BENT.bio  = b_mean(:,end);

save([fpath 'Hist_ctrlclim_All_fishobs_',mod,'_Last_mo_' cfile '.mat'],'Sml_f','Sml_p','Sml_d',... 
    'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')
save([fpath 'Hist_obsclim_All_fishobs_',mod,'_Last_mo_' cfile '.mat'],'Sml_f','Sml_p','Sml_d',... 
    'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')


