% FEISTY output at all locations

clear
close all

%/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/.../OneDeg
%Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100

%/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/.../OneDeg
%Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];

mods = {'assessment','effective','nominal','FFmsy_creep','FFmsy_nominal','FFmsymax_creep','FFmsymax_nominal',...
    'FFmsymin_creep','FFmsymin_nominal'};

%%
for i=4:length(mods)

    mod = ['All_fish_obs_' mods{i}];
    mod2 = ['All_fishobs_' mods{i}];

    %% SP
    ncid = netcdf.open([fpath 'Spinup_ctrlclim_',mod,'_empHP_sml_p.nc'],'NC_NOWRITE');
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
    ncid = netcdf.open([fpath 'Spinup_ctrlclim_',mod,'_empHP_sml_f.nc'],'NC_NOWRITE');
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
    ncid = netcdf.open([fpath 'Spinup_ctrlclim_',mod,'_empHP_sml_d.nc'],'NC_NOWRITE');
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
    ncid = netcdf.open([fpath 'Spinup_ctrlclim_',mod,'_empHP_med_p.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MP.bio = biomass;
    MP.yield = yield;
    clear yield
    clear biomass

    % MF
    ncid = netcdf.open([fpath 'Spinup_ctrlclim_',mod,'_empHP_med_f.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MF.bio = biomass;
    MF.yield = yield;
    clear yield
    clear biomass

    % MD
    ncid = netcdf.open([fpath 'Spinup_ctrlclim_',mod,'_empHP_med_d.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MD.bio = biomass;
    MD.yield = yield;
    clear yield
    clear biomass

    % LP
    ncid = netcdf.open([fpath 'Spinup_ctrlclim_',mod,'_empHP_lrg_p.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    LP.bio = biomass;
    LP.yield = yield;
    clear yield
    clear biomass

    % LD
    ncid = netcdf.open([fpath 'Spinup_ctrlclim_',mod,'_empHP_lrg_d.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    LD.bio = biomass;
    LD.yield = yield;
    clear yield
    clear biomass

    % Benthic material
    ncid = netcdf.open([fpath 'Spinup_ctrlclim_',mod,'_empHP_bent.nc'],'NC_NOWRITE');
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
    sp_tmean=mean(SP.bio,1,'omitnan');
    sf_tmean=mean(SF.bio,1,'omitnan');
    sd_tmean=mean(SD.bio,1,'omitnan');
    mp_tmean=mean(MP.bio,1,'omitnan');
    mf_tmean=mean(MF.bio,1,'omitnan');
    md_tmean=mean(MD.bio,1,'omitnan');
    lp_tmean=mean(LP.bio,1,'omitnan');
    ld_tmean=mean(LD.bio,1,'omitnan');
    b_tmean=mean(Bent.bio,1,'omitnan');

    mp_tmcatch=mean(MP.yield,1,'omitnan');
    mf_tmcatch=mean(MF.yield,1,'omitnan');
    md_tmcatch=mean(MD.yield,1,'omitnan');
    lp_tmcatch=mean(LP.yield,1,'omitnan');
    ld_tmcatch=mean(LD.yield,1,'omitnan');

    %% Space
    sp_mean=mean(SP.bio(:,nt-11:nt),2,'omitnan');
    sf_mean=mean(SF.bio(:,nt-11:nt),2,'omitnan');
    sd_mean=mean(SD.bio(:,nt-11:nt),2,'omitnan');
    mp_mean=mean(MP.bio(:,nt-11:nt),2,'omitnan');
    mf_mean=mean(MF.bio(:,nt-11:nt),2,'omitnan');
    md_mean=mean(MD.bio(:,nt-11:nt),2,'omitnan');
    lp_mean=mean(LP.bio(:,nt-11:nt),2,'omitnan');
    ld_mean=mean(LD.bio(:,nt-11:nt),2,'omitnan');
    b_mean =mean(Bent.bio(:,nt-11:nt),2,'omitnan');

    mp_mcatch=mean(MP.yield(:,nt-11:nt),2,'omitnan');
    mf_mcatch=mean(MF.yield(:,nt-11:nt),2,'omitnan');
    md_mcatch=mean(MD.yield(:,nt-11:nt),2,'omitnan');
    lp_mcatch=mean(LP.yield(:,nt-11:nt),2,'omitnan');
    ld_mcatch=mean(LD.yield(:,nt-11:nt),2,'omitnan');

    %%
    save([fpath 'Means_Spinup_ctrlclim_',mod,'_' cfile '.mat'],'time',...
        'sf_tmean','sp_tmean','sd_tmean',...
        'mf_tmean','mp_tmean','md_tmean',...
        'lp_tmean','ld_tmean','b_tmean',...
        'mf_tmcatch','mp_tmcatch','md_tmcatch',...
        'lp_tmcatch','ld_tmcatch',...
        'sf_mean','sp_mean','sd_mean',...
        'mf_mean','mp_mean','md_mean',...
        'lp_mean','ld_mean','b_mean',...
        'mf_mcatch','mp_mcatch','md_mcatch',...
        'lp_mcatch','ld_mcatch')


    % Save last year for initializing forecast runs
    Sml_f.bio = mean(SF.bio(:,nt-11:nt),2,'omitnan');
    Sml_p.bio = mean(SP.bio(:,nt-11:nt),2,'omitnan');
    Sml_d.bio = mean(SD.bio(:,nt-11:nt),2,'omitnan');
    Med_f.bio = mean(MF.bio(:,nt-11:nt),2,'omitnan');
    Med_p.bio = mean(MP.bio(:,nt-11:nt),2,'omitnan');
    Med_d.bio = mean(MD.bio(:,nt-11:nt),2,'omitnan');
    Lrg_p.bio = mean(LP.bio(:,nt-11:nt),2,'omitnan');
    Lrg_d.bio = mean(LD.bio(:,nt-11:nt),2,'omitnan');
    BENT.bio  = mean(Bent.bio(:,nt-11:nt),2,'omitnan');

    save([fpath 'PI_ctrlclim_',mod,'_Last_mo_' cfile '.mat'],...
        'Sml_f','Sml_p','Sml_d',...
        'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')
    save([fpath 'PI_ctrlclim_',mod2,'_Last_mo_' cfile '.mat'],...
        'Sml_f','Sml_p','Sml_d',...
        'Med_f','Med_p','Med_d','Lrg_p','Lrg_d','BENT')

end
