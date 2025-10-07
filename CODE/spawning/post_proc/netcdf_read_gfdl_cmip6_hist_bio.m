% FEISTY output at all locations

clear 
close all

%%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

%fpath=['/project/Feisty/NC/spawning/' cfile '/CMIP6/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/spawning/' cfile '/CMIP6/'];

mod = 'Historic_const_spawning_All_fish03_1950_';

%% SP
ncid = netcdf.open([fpath mod 'sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%%
[ni,nt] = size(biomass);

SP.bio = biomass;
clear biomass 

%% SF
ncid = netcdf.open([fpath mod 'sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
clear biomass 

% SD
ncid = netcdf.open([fpath mod 'sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
clear biomass

% MP
ncid = netcdf.open([fpath mod 'med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
clear biomass 

% MF
ncid = netcdf.open([fpath mod 'med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
clear biomass 

% MD
ncid = netcdf.open([fpath mod 'med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
clear biomass 

% LP
ncid = netcdf.open([fpath mod 'lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
clear biomass 

% LD
ncid = netcdf.open([fpath mod 'lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
clear biomass 

% Benthic material
ncid = netcdf.open([fpath mod 'bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
clear biomass


%% Take means for visualization

%Time
sp_tmean = mean(SP.bio,1,'omitnan');
sf_tmean = mean(SF.bio,1,'omitnan');
sd_tmean = mean(SD.bio,1,'omitnan');
mp_tmean = mean(MP.bio,1,'omitnan');
mf_tmean = mean(MF.bio,1,'omitnan');
md_tmean = mean(MD.bio,1,'omitnan');
lp_tmean = mean(LP.bio,1,'omitnan');
ld_tmean = mean(LD.bio,1,'omitnan');
b_tmean  = mean(Bent.bio,1,'omitnan');

%% Space
sp_sbio = mean(SP.bio,2,'omitnan');
sf_sbio = mean(SF.bio,2,'omitnan');
sd_sbio = mean(SD.bio,2,'omitnan');
mp_sbio = mean(MP.bio,2,'omitnan');
mf_sbio = mean(MF.bio,2,'omitnan');
md_sbio = mean(MD.bio,2,'omitnan');
lp_sbio = mean(LP.bio,2,'omitnan');
ld_sbio = mean(LD.bio,2,'omitnan');
b_sbio  = mean(Bent.bio,2,'omitnan');

%% Annual means
nyr = nt/12;
% st=1:12:length(time);
% en=12:12:length(time);
st=1:12:nt;
en=12:12:nt;

for n=1:length(st)
    % mean biomass
    sp_abio(:,n)=mean(SP.bio(:,st(n):en(n)),2,'omitnan');
    sf_abio(:,n)=mean(SF.bio(:,st(n):en(n)),2,'omitnan');
    sd_abio(:,n)=mean(SD.bio(:,st(n):en(n)),2,'omitnan');
    mp_abio(:,n)=mean(MP.bio(:,st(n):en(n)),2,'omitnan');
    mf_abio(:,n)=mean(MF.bio(:,st(n):en(n)),2,'omitnan');
    md_abio(:,n)=mean(MD.bio(:,st(n):en(n)),2,'omitnan');
    lp_abio(:,n)=mean(LP.bio(:,st(n):en(n)),2,'omitnan');
    ld_abio(:,n)=mean(LD.bio(:,st(n):en(n)),2,'omitnan');
    b_abio(:,n)=mean(Bent.bio(:,st(n):en(n)),2,'omitnan');

end

%%
save([fpath 'Time_Means_' mod cfile '.mat'],'time',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean')

save([fpath 'Space_Means_' mod cfile '.mat'],'time',...
    'sf_sbio','sp_sbio','sd_sbio',...
    'mf_sbio','mp_sbio','md_sbio',...
    'lp_sbio','ld_sbio','b_sbio')

save([fpath 'Annual_Means_' mod cfile '.mat'],'time',...
    'sf_abio','sp_abio','sd_abio',...
    'mf_abio','mp_abio','md_abio',...
    'lp_abio','ld_abio','b_abio')

%%
% mo = time/12;
mo = time;

figure
plot(mo,log10(lp_tmean),'b'); hold on;
plot(mo,log10(mf_tmean),'r'); hold on;
plot(mo,log10(ld_tmean),'k'); hold on;

