% FEISTY output at all locations

clear
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

esms = {'IPSL','UKESM','CESM2-WACCM'};

for m=2:3 %1:length(esms)

    mod = esms{m};

    exper = [mod '_ssp534-over_pristine'];
    
    fpath=['/project/Feisty/NC/WG2300/',cfile,'/',mod,'/'];

    %% SP
    ncid = netcdf.open([fpath exper '_empHP_sml_p.nc'],'NC_NOWRITE');
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
    ncid = netcdf.open([fpath exper '_empHP_sml_f.nc'],'NC_NOWRITE');
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
    ncid = netcdf.open([fpath exper '_empHP_sml_d.nc'],'NC_NOWRITE');
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
    ncid = netcdf.open([fpath exper '_empHP_med_p.nc'],'NC_NOWRITE');
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
    ncid = netcdf.open([fpath exper '_empHP_med_f.nc'],'NC_NOWRITE');
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
    ncid = netcdf.open([fpath exper '_empHP_med_d.nc'],'NC_NOWRITE');
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
    ncid = netcdf.open([fpath exper '_empHP_lrg_p.nc'],'NC_NOWRITE');
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
    ncid = netcdf.open([fpath exper '_empHP_lrg_d.nc'],'NC_NOWRITE');
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
    ncid = netcdf.open([fpath exper '_empHP_bent.nc'],'NC_NOWRITE');
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

    %Space
    t=time;
    mo=t/12;
    mo=mo+2040;
    yr1=find(mo>2099 & mo<=2100);
    yr3=find(mo>2299 & mo<=2300);

    sp_mean1=mean(SP.bio(:,yr1),2,'omitnan');
    sf_mean1=mean(SF.bio(:,yr1),2,'omitnan');
    sd_mean1=mean(SD.bio(:,yr1),2,'omitnan');
    mp_mean1=mean(MP.bio(:,yr1),2,'omitnan');
    mf_mean1=mean(MF.bio(:,yr1),2,'omitnan');
    md_mean1=mean(MD.bio(:,yr1),2,'omitnan');
    lp_mean1=mean(LP.bio(:,yr1),2,'omitnan');
    ld_mean1=mean(LD.bio(:,yr1),2,'omitnan');
    b_mean1 =mean(Bent.bio(:,yr1),2,'omitnan');

    sp_mean3=mean(SP.bio(:,yr3),2,'omitnan');
    sf_mean3=mean(SF.bio(:,yr3),2,'omitnan');
    sd_mean3=mean(SD.bio(:,yr3),2,'omitnan');
    mp_mean3=mean(MP.bio(:,yr3),2,'omitnan');
    mf_mean3=mean(MF.bio(:,yr3),2,'omitnan');
    md_mean3=mean(MD.bio(:,yr3),2,'omitnan');
    lp_mean3=mean(LP.bio(:,yr3),2,'omitnan');
    ld_mean3=mean(LD.bio(:,yr3),2,'omitnan');
    b_mean3 =mean(Bent.bio(:,yr3),2,'omitnan');

    if m==3
        exper2 = [mod '_ssp534_zooc_pristine'];
    else
        exper2 = [mod '_ssp534_pristine'];
    end

    save([fpath 'Means_' exper2 '_' cfile '.mat'],'time','mo',...
        'sf_tmean','sp_tmean','sd_tmean',...
        'mf_tmean','mp_tmean','md_tmean',...
        'lp_tmean','ld_tmean','b_tmean',...
        'sf_mean1','sp_mean1','sd_mean1',...
        'mf_mean1','mp_mean1','md_mean1',...
        'lp_mean1','ld_mean1','b_mean1',...
        'sf_mean3','sp_mean3','sd_mean3',...
        'mf_mean3','mp_mean3','md_mean3',...
        'lp_mean3','ld_mean3','b_mean3')


    %%
    % figure
    % plot(time,mf_tmean,'r'); hold on
    % plot(time,lp_tmean,'b'); hold on
    % plot(time,ld_tmean,'k')

    %% Fish-MIP OUTPUTS =================================================

    % PREFERRED (all units = gWW/m2)
    
    %total pelagic (Linf <30cm) biomass bp30cm = 360x180xMOs
    SPel = SF.bio + MF.bio;

    %total pelagic (>=30 cm and <90cm) biomass bp30to90cm = 360x180xMOs
    %none

    %total pelagic (>=90cm) biomass bp90cm = 360x180xMOs
    LPel = SP.bio + MP.bio + LP.bio;

    %total pelagic biomass tpb = 360x180xMOs
    allPel = SPel + LPel;

    %total demersal biomass tdb = 360x180xMOs
    allD = SD.bio + MD.bio + LD.bio;

    %total consumber biomass tcb = 360x180xMOs
    % vtcb =
    allC = allPel + allD + Bent.bio;

    %total consumber biomass in log10 bins tcblog10 = 360x180xMOsx6
    %(1g, 10g, 100g, 1kg, 10kg, 100kg)
    %Small <1g      (0.001-0.5g; 0.02g mean)    (0.46-3.68cm; mean 1.3cm)
    %Med 1-100g     (0.5-250g; 11.2g mean)      (3.68-29.24cm; mean 10.4cm)
    %Lrg 1-100kg    (250-125000g; 5600g mean)   (29.24-232.08cm; mean 82.4cm)

    % SECONDARY

    %total demersal (Linf <30cm) biomass bd30cm = 360x180xMOs
    %none

    %total demersal (>=30 cm and <90cm) biomass bd30to90cm = 360x180xMOs
    %none

    %total demersal (>=90cm) biomass bd90cm = 360x180xMOs
    %LDem = allD;
 
    save([fpath 'FishMIP_outputs_monthly_' exper2 '_' cfile '.mat'],'time','mo',...
        'allPel','allD','allC','SPel','LPel');

end