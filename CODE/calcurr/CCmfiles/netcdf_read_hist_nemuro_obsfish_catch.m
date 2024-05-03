% FEISTY output at all locations

clear
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';
esm1 = {'IPSL','GFDL','HAD'};
esm2 = {'ipsl','gfdl','hadley'};
harv = 'All_fishobs';

%%
for z=1:2

    vers = esm1{z};
    vers2 = esm2{z};

    fpath=['/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/',vers,'/'];
    Cdir = ['/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/',vers,'down/'];

    load([Cdir 'Data_grid_nemuro_',vers2,'.mat'],'GRD')

    %% Time from Benthos
    ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_bent.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    nt = length(time);

    clear biomass

    % MP
    ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_med_p.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MP.yield = yield;
    clear biomass yield

    % MF
    ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_med_f.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MF.yield = yield;
    clear biomass yield

    % MD
    ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_med_d.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MD.yield = yield;
    clear biomass yield

    % LP
    ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_lrg_p.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    LP.yield = yield;
    clear biomass yield

    % LD
    ncid = netcdf.open([fpath 'Hist_' vers '_' harv '_lrg_d.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    LD.yield = yield;
    clear biomass yield

    %% Total yields = area sums
    % Area in m2 instead of km2, matrix
    aream2 = GRD.AREA .* 1e6;
    aream2_mat = repmat(aream2,1,nt);

    mf_ay= (MF.yield .* aream2_mat);
    mp_ay= (MP.yield .* aream2_mat);
    md_ay= (MD.yield .* aream2_mat);
    lp_ay= (LP.yield .* aream2_mat);
    ld_ay= (LD.yield .* aream2_mat);

    %% Time total
    mf_tty=sum(mf_ay,1,"omitnan");
    mp_tty=sum(mp_ay,1,"omitnan");
    md_tty=sum(md_ay,1,"omitnan");
    lp_tty=sum(lp_ay,1,"omitnan");
    ld_tty=sum(ld_ay,1,"omitnan");

    %% Space total
    mf_sty=sum(mf_ay,2,"omitnan");
    mp_sty=sum(mp_ay,2,"omitnan");
    md_sty=sum(md_ay,2,"omitnan");
    lp_sty=sum(lp_ay,2,"omitnan");
    ld_sty=sum(ld_ay,2,"omitnan");

    %% Each year
    a = 1:12:nt; % start of each yr
    b = 12:12:nt; % end of each yr
    ny = length(a);
    nid = length(GRD.ID);

    mB = NaN*ones(nid,ny);
    ytcMF = mB;
    ytcMP = mB;
    ytcMD = mB;
    ytcLP = mB;
    ytcLD = mB;
    for i = 1:(nt/12)
        ytcMF(:,i) = sum(mf_ay(:,a(i):b(i)),2,"omitnan");
        ytcMP(:,i) = sum(mf_ay(:,a(i):b(i)),2,"omitnan");
        ytcMD(:,i) = sum(mf_ay(:,a(i):b(i)),2,"omitnan");
        ytcLP(:,i) = sum(mf_ay(:,a(i):b(i)),2,"omitnan");
        ytcLD(:,i) = sum(mf_ay(:,a(i):b(i)),2,"omitnan");
    end

    %%
    save([fpath 'Means_Hist_' vers '_' harv '_catch_' cfile '.mat'],...
        'mf_tty','mp_tty','md_tty','lp_tty','ld_tty',...
        'mf_sty','mp_sty','md_sty','lp_sty','ld_sty',...
        'ytcMF','ytcMP','ytcMD','ytcLP','ytcLD');


end

