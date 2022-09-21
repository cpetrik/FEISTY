% FEISTY output at all locations

clear all
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];

%% intermediate cycles
% load([fpath 'Spinup_ctrlclim_pristine_end_cycle_1.mat'])
%
% % Save last year for initializing forecast runs
% Sml_f.bio = ySF(:,1);
% Sml_p.bio = ySP(:,1);
% Sml_d.bio = ySD(:,1);
% Med_f.bio = yMF(:,1);
% Med_p.bio = yMP(:,1);
% Med_d.bio = yMD(:,1);
% Lrg_p.bio = yLP(:,1);
% Lrg_d.bio = yLD(:,1);
% BENT.mass = yB(:,1);
%
% save([fpath 'Last_mo_spinup_ctrlclim_pristine_' cfile '.mat'],'Sml_f','Sml_p','Sml_d',...
%     'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')

%% loop over cycles
nid = 670589;
ntc = 240;
nta = ntc*6;
yrs = 1841:1960;
st = 1:ntc:nta;
en = ntc:ntc:nta;

SP.bio = nan*ones(nid,nta);

%%
for c=1:10

    %% SP
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_sml_p_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    %%
    SP.bio(:,st(c):en(c)) = biomass;
    clear biomass

    %% SF
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_sml_f_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    SF.bio(:,st(c):en(c)) = biomass;
    clear biomass

    % SD
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_sml_d_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    SD.bio(:,st(c):en(c)) = biomass;
    clear biomass

    % MP
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_med_p_cycle',num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);...

    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MP.bio(:,st(c):en(c)) = biomass;
    clear biomass

    % MF
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_med_f_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MF.bio(:,st(c):en(c)) = biomass;
    clear biomass

    % MD
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_med_d_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MD.bio(:,st(c):en(c)) = biomass;
    clear biomass

    % LP
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_lrg_p_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    LP.bio(:,st(c):en(c)) = biomass;
    clear biomass

    % LD
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_lrg_d_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    LD.bio(:,st(c):en(c)) = biomass;
    clear biomass

    % Benthic material
    ncid = netcdf.open([fpath 'PI_ctrlclim_pristine_empHP_bent_cycle',...
        num2str(c),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    Bent.bio(:,st(c):en(c)) = biomass;
    clear biomass

end

%% Take means

%Time
sp_tmean=mean(SP.bio,1,'omitnan');
sf_tmean=mean(SF.bio,1,'omitnan');
sd_tmean=mean(SD.bio,1,'omitnan');
mp_tmean=mean(MP.bio,1,'omitnan');
mf_tmean=mean(MF.bio,1,'omitnan');
md_tmean=mean(MD.bio,1,'omitnan');
lp_tmean=mean(LP.bio,1,'omitnan');
ld_tmean=mean(LD.bio,1,'omitnan');
b_tmean=mean(Bent.bio,1,'omitnan');

%Space
sp_mean=smean(SP.bio,2,'omitnan');
sf_mean=smean(SF.bio,2,'omitnan');
sd_mean=smean(SD.bio,2,'omitnan');
mp_mean=smean(MP.bio,2,'omitnan');
mf_mean=smean(MF.bio,2,'omitnan');
md_mean=smean(MD.bio,2,'omitnan');
lp_mean=smean(LP.bio,2,'omitnan');
ld_mean=smean(LD.bio,2,'omitnan');
b_mean =smean(Bent.bio,2,'omitnan');

%% Each year
st=1:12:length(time);
en=12:12:length(time);

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

save([fpath 'Means_PI_ctrlclim_pristine_',cfile,'.mat'],'time',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean',...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean',...
    'sf_smean','sp_smean','sd_smean',...
    'mf_smean','mp_smean','md_smean',...
    'lp_smean','ld_smean','b_smean')

%% Save last year for initializing forecast run
% Sml_f.bio = mean(SF.bio(:,nta-11:nta),2,'omitnan');
% Sml_p.bio = mean(SP.bio(:,nta-11:nta),2,'omitnan');
% Sml_d.bio = mean(SD.bio(:,nta-11:nta),2,'omitnan');
% Med_f.bio = mean(MF.bio(:,nta-11:nta),2,'omitnan');
% Med_p.bio = mean(MP.bio(:,nta-11:nta),2,'omitnan');
% Med_d.bio = mean(MD.bio(:,nta-11:nta),2,'omitnan');
% Lrg_p.bio = mean(LP.bio(:,nta-11:nta),2,'omitnan');
% Lrg_d.bio = mean(LD.bio(:,nta-11:nta),2,'omitnan');
% BENT.bio  = mean(Bent.bio(:,nta-11:nta),2,'omitnan');
% 
% save([fpath 'Hist_ctrlclim_pristine_Last_mo_cycle',num2str(c),'_',cfile '.mat'],'Sml_f','Sml_p','Sml_d',...
%     'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')
% save([fpath 'Hist_obsclim_pristine_Last_mo_cycle',num2str(c),'_', cfile '.mat'],'Sml_f','Sml_p','Sml_d',...
%     'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')





