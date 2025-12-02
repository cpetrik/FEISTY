% CORE-forced run with movement
% used IC of Spinup for 20 yrs
% catch output 1988-2007

clear
close all

%%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_CC80_RE00100';
%fpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/CORE/'];
fpath=['/project/Feisty/NC/Matlab_new_size/' cfile '/CORE/'];

mods = {'enc','ingest','mort','nu','preyconc'};

for m =1:length(mods)

    exper = ['CORE_Hindcast_move_',mods{m},'_v28_dt12h_All_fish03_spin20IC_yield' ];

    %% MP
    ncid = netcdf.open([fpath exper '_med_p.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MP.yield = yield;

    clear yield

    % MF
    ncid = netcdf.open([fpath exper '_med_f.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MF.yield = yield;

    clear yield

    % MD
    ncid = netcdf.open([fpath exper '_med_d.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MD.yield = yield;

    clear yield

    % LP
    ncid = netcdf.open([fpath exper '_lrg_p.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    LP.yield = yield;

    clear yield

    % LD
    ncid = netcdf.open([fpath exper '_lrg_d.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    LD.yield = yield;

    clear  yield


    %% Take means and totals
    % Catch totals gridded
    vpath = '/project/Feisty/GCM_Data/CORE-forced/';

    %1-D
    load([vpath 'Data_grid_ocean_cobalt_ESM2Mcore.mat'],'GRD');
    GRD1 = GRD;
    clear GRD

    %2-D
    load([vpath 'Data_hindcast_grid_cp2D.mat'],'GRD')
    GRD2 = GRD;
    clear GRD

    ID = GRD1.ID;

    %AREA units 'm^2'
    [ni,nj]  = size(GRD2.area);
    [nid,nt] = size(LD.yield);

    MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
    nyr = nt/12;
    mos = repmat(MNTH,nid,nyr);
    mns = repmat(MNTH,nid,1);
    area_mat = repmat(GRD1.area,1,nt);

    %% Units
    units_yield = 'g_m2_day';
    units_catch = 'g_mo';

    %                  #days/mo   total m2
    MF.catch = MF.yield .*mos .*area_mat;
    MP.catch = MP.yield .*mos .*area_mat;
    MD.catch = MD.yield .*mos .*area_mat;
    LP.catch = LP.yield .*mos .*area_mat;
    LD.catch = LD.yield .*mos .*area_mat;

    %% Time
    %mean yield per mo
    mf_tmy=mean(MF.yield,1,'omitnan');
    mp_tmy=mean(MP.yield,1,'omitnan');
    md_tmy=mean(MD.yield,1,'omitnan');
    lp_tmy=mean(LP.yield,1,'omitnan');
    ld_tmy=mean(LD.yield,1,'omitnan');

    %mean catch per mo
    mf_tmc=mean(MF.catch,1,'omitnan');
    mp_tmc=mean(MP.catch,1,'omitnan');
    md_tmc=mean(MD.catch,1,'omitnan');
    lp_tmc=mean(LP.catch,1,'omitnan');
    ld_tmc=mean(LD.catch,1,'omitnan');

    %total yield per mo
    mf_tty=sum(MF.yield,1,'omitnan');
    mp_tty=sum(MP.yield,1,'omitnan');
    md_tty=sum(MD.yield,1,'omitnan');
    lp_tty=sum(LP.yield,1,'omitnan');
    ld_tty=sum(LD.yield,1,'omitnan');

    %total catch per mo
    mf_ttc=sum(MF.catch,1,'omitnan');
    mp_ttc=sum(MP.catch,1,'omitnan');
    md_ttc=sum(MD.catch,1,'omitnan');
    lp_ttc=sum(LP.catch,1,'omitnan');
    ld_ttc=sum(LD.catch,1,'omitnan');


    %% Spatially
    %mean yield per mo
    mf_smy=mean(MF.yield,2,'omitnan');
    mp_smy=mean(MP.yield,2,'omitnan');
    md_smy=mean(MD.yield,2,'omitnan');
    lp_smy=mean(LP.yield,2,'omitnan');
    ld_smy=mean(LD.yield,2,'omitnan');

    %mean catch per mo
    mf_smc=mean(MF.catch,2,'omitnan');
    mp_smc=mean(MP.catch,2,'omitnan');
    md_smc=mean(MD.catch,2,'omitnan');
    lp_smc=mean(LP.catch,2,'omitnan');
    ld_smc=mean(LD.catch,2,'omitnan');

    %total yield per mo
    mf_sty=sum(MF.yield,2,'omitnan');
    mp_sty=sum(MP.yield,2,'omitnan');
    md_sty=sum(MD.yield,2,'omitnan');
    lp_sty=sum(LP.yield,2,'omitnan');
    ld_sty=sum(LD.yield,2,'omitnan');

    %total catch per mo
    mf_stc=sum(MF.catch,2,'omitnan');
    mp_stc=sum(MP.catch,2,'omitnan');
    md_stc=sum(MD.catch,2,'omitnan');
    lp_stc=sum(LP.catch,2,'omitnan');
    ld_stc=sum(LD.catch,2,'omitnan');


    %% Every year
    st=1:12:nt;
    en=12:12:nt;

    mf_tac = nan*ones(nid,nyr);
    mp_tac = nan*ones(nid,nyr);
    md_tac = nan*ones(nid,nyr);
    lp_tac = nan*ones(nid,nyr);
    ld_tac = nan*ones(nid,nyr);
    mf_may = nan*ones(nid,nyr);
    mp_may = nan*ones(nid,nyr);
    md_may = nan*ones(nid,nyr);
    lp_may = nan*ones(nid,nyr);
    ld_may = nan*ones(nid,nyr);

    Cmf=NaN*ones(ni,nj,nyr);
    Cmp=NaN*ones(ni,nj,nyr);
    Cmd=NaN*ones(ni,nj,nyr);
    Clp=NaN*ones(ni,nj,nyr);
    Cld=NaN*ones(ni,nj,nyr);


    for n=1:length(st)

        mp_tac(:,n)=sum(MP.catch(:,st(n):en(n)),2,'omitnan');
        mf_tac(:,n)=sum(MF.catch(:,st(n):en(n)),2,'omitnan');
        md_tac(:,n)=sum(MD.catch(:,st(n):en(n)),2,'omitnan');
        lp_tac(:,n)=sum(LP.catch(:,st(n):en(n)),2,'omitnan');
        ld_tac(:,n)=sum(LD.catch(:,st(n):en(n)),2,'omitnan');

        mp_may(:,n)=mean(MP.yield(:,st(n):en(n)),2,'omitnan');
        mf_may(:,n)=mean(MF.yield(:,st(n):en(n)),2,'omitnan');
        md_may(:,n)=mean(MD.yield(:,st(n):en(n)),2,'omitnan');
        lp_may(:,n)=mean(LP.yield(:,st(n):en(n)),2,'omitnan');
        ld_may(:,n)=mean(LD.yield(:,st(n):en(n)),2,'omitnan');

        Zmf=NaN*ones(ni,nj);
        Zmp=NaN*ones(ni,nj);
        Zmd=NaN*ones(ni,nj);
        Zlp=NaN*ones(ni,nj);
        Zld=NaN*ones(ni,nj);

        Zmf(ID)=squeeze(mf_tac(:,n));
        Zmp(ID)=squeeze(mp_tac(:,n));
        Zmd(ID)=squeeze(md_tac(:,n));
        Zlp(ID)=squeeze(lp_tac(:,n));
        Zld(ID)=squeeze(ld_tac(:,n));

        Cmf(:,:,n) = Zmf;
        Cmp(:,:,n) = Zmp;
        Cmd(:,:,n) = Zmd;
        Clp(:,:,n) = Zlp;
        Cld(:,:,n) = Zld;
    end

    tmn = mf_tac + mp_tac + md_tac + lp_tac + ld_tac;
    stmn = sum(tmn,'omitnan');

    mp_tsac = sum(mp_tac,'omitnan');
    mf_tsac = sum(mf_tac,'omitnan');
    md_tsac = sum(md_tac,'omitnan');
    lp_tsac = sum(lp_tac,'omitnan');
    ld_tsac = sum(ld_tac,'omitnan');

    units_tac = 'g_yr';


    %%
    save([fpath 'Means_' exper '_' cfile '.mat'],'units_yield','units_catch',...
        'mf_tmy','mp_tmy','md_tmy','lp_tmy','ld_tmy',...
        'mf_tmc','mp_tmc','md_tmc','lp_tmc','ld_tmc',...
        'mf_tty','mp_tty','md_tty','lp_tty','ld_tty',...
        'mf_ttc','mp_ttc','md_ttc','lp_ttc','ld_ttc',...
        'mf_smy','mp_smy','md_smy','lp_smy','ld_smy',...
        'mf_smc','mp_smc','md_smc','lp_smc','ld_smc',...
        'mf_sty','mp_sty','md_sty','lp_sty','ld_sty',...
        'mf_stc','mp_stc','md_stc','lp_stc','ld_stc',...
        'units_tac');

    save([fpath 'Annual_Means_' exper '_' cfile '.mat'],...
        'mf_tac','mp_tac','md_tac','lp_tac','ld_tac',...
        'mf_tsac','mp_tsac','md_tsac','lp_tsac','ld_tsac',...
        'units_yield','units_catch','units_tac',...
        'mf_may','mp_may','md_may','lp_may','ld_may',...
        'Cmf','Cmp','Cmd','Clp','Cld');


end
