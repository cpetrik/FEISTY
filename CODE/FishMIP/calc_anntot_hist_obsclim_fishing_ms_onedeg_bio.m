% FEISTY output at all locations
% Years 1960-2010 will be good
% annual means for biomass and annual sums for catches are okay

clear
close all

%%
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];

mods = {'All_fishobs_assessment','All_fishobs_effective','All_fishobs_nominal',...
    'All_fish_obs_FFmsy_creep','All_fish_obs_FFmsy_nominal','All_fish_obs_FFmsymax_creep',...
    'All_fish_obs_FFmsymax_nominal','All_fish_obs_FFmsymin_creep','All_fish_obs_FFmsymin_nominal'};

%% Area for totals
cpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';
load([cpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat']);
load([cpath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat'], 'GRD');
load([cpath 'cellarea_onedeg.mat']);

ID = GRD.ID;

%area units 'm^2'
area = cell_area(ID);

%%
for i=4:length(mods)
    alt = mods{i};

    %% MF
    ncid = netcdf.open([fpath 'Hist_obsclim_',alt,'_empHP_med_f.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MF.yield = yield;
    clear yield
    clear biomass

    %% MP
    ncid = netcdf.open([fpath 'Hist_obsclim_',alt,'_empHP_med_p.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MP.yield = yield;
    clear yield
    clear biomass

    % MD
    ncid = netcdf.open([fpath 'Hist_obsclim_',alt,'_empHP_med_d.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MD.yield = yield;
    clear yield
    clear biomass

    % LP
    ncid = netcdf.open([fpath 'Hist_obsclim_',alt,'_empHP_lrg_p.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    LP.yield = yield;
    clear yield
    clear biomass

    % LD
    ncid = netcdf.open([fpath 'Hist_obsclim_',alt,'_empHP_lrg_d.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    LD.yield = yield;
    clear yield
    clear biomass

    %% Units

    [nid,nt] = size(MF.yield);

    MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
    nyr = nt/12;
    mos = repmat(MNTH,nid,nyr);
    mns = repmat(MNTH,nid,1);
    area_mat = repmat(area,1,nt);

    units_yield = 'g_m2_day';
    units_catch = 'g_mo';

    MF.catch = MF.yield .*mos .*area_mat;
    MP.catch = MP.yield .*mos .*area_mat;
    MD.catch = MD.yield .*mos .*area_mat;
    LP.catch = LP.yield .*mos .*area_mat;
    LD.catch = LD.yield .*mos .*area_mat;

    %% Time total catch per mo
    mf_ttc=sum(MF.catch,1,'omitnan');
    mp_ttc=sum(MP.catch,1,'omitnan');
    md_ttc=sum(MD.catch,1,'omitnan');
    lp_ttc=sum(LP.catch,1,'omitnan');
    ld_ttc=sum(LD.catch,1,'omitnan');

    %% Spatially total catch per mo
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
    for n=1:length(st)

        mp_tac(:,n)=sum(MP.catch(:,st(n):en(n)),2,'omitnan');
        mf_tac(:,n)=sum(MF.catch(:,st(n):en(n)),2,'omitnan');
        md_tac(:,n)=sum(MD.catch(:,st(n):en(n)),2,'omitnan');
        lp_tac(:,n)=sum(LP.catch(:,st(n):en(n)),2,'omitnan');
        ld_tac(:,n)=sum(LD.catch(:,st(n):en(n)),2,'omitnan');
    end

    tmn = mf_tac + mp_tac + md_tac + lp_tac + ld_tac;
    stmn = sum(tmn);

    mp_tsac = sum(mp_tac);
    mf_tsac = sum(mf_tac);
    md_tsac = sum(md_tac);
    lp_tsac = sum(lp_tac);
    ld_tsac = sum(ld_tac);

    %% Quick check
    % 1 tonne = 1 MT = 1e-6 g
    figure
    plot(1:50,mf_tsac*1e-6,'r'); hold on;
    plot(1:50,lp_tsac*1e-6,'b'); hold on;
    plot(1:50,ld_tsac*1e-6,'k'); hold on;

    %% Save
    save([fpath 'Catch_anntots_Hist_obsclim_',alt,'_',cfile,'.mat'],'units_catch',...
        'mf_ttc','mp_ttc','md_ttc',...
        'lp_ttc','ld_ttc',...
        'mf_tac','mp_tac','md_tac',...
        'lp_tac','ld_tac',...
        'mf_stc','mp_stc','md_stc',...
        'lp_stc','ld_stc');

end



