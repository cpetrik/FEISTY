% 0.5 degree global MOM6-COBALTv3-FEISTY only


clear
close all

%fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';
fpath = '/project/Feisty/NC/Global_COBALT_FEISTY/';

%gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
gpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%%
%load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet','z_l_units','z_l_long_name','z_l')

%%
ncdisp([fpath 'ocean_cobalt_tracers_month_z.199001-199412.ndet.nc'])

%%
st = 1990:5:2019;
en = 1994:5:2019;

tSP = nan*ones(1,60*length(st));
tLP = tSP;
tSZ = tSP;
tDE = tSP;
tDI = tSP;

sSP = nan*ones(720,576,length(st));
sLP = sSP;
sSZ = sSP;
sDE = sSP;
sDI = sSP;

vSP = nan*ones(35,length(st));
vLP = vSP;
vSZ = vSP;
vDE = vSP;
vDI = vSP;

%%
for y=1%:6

    %% ndet
    ncid = netcdf.open([fpath 'ocean_cobalt_tracers_month_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.ndet.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 5
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    ndet(ndet>1e19) = nan;

    %% nsm
    ncid = netcdf.open([fpath 'ocean_cobalt_tracers_month_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.nsm.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 5
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    nsm(nsm>1e19) = nan;

    %% nlg
    ncid = netcdf.open([fpath 'ocean_cobalt_tracers_month_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.nlg.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 5
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    nlg(nlg>1e19) = nan;

    %% ndiaz
    ncid = netcdf.open([fpath 'ocean_cobalt_tracers_month_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.ndi.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 5
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    ndi(ndi>1e19) = nan;

    %% SZ
    ncid = netcdf.open([fpath 'ocean_cobalt_tracers_month_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.nsmz.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 5
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    nsmz(nsmz>1e19) = nan;

    %% thkcello
    load([gpath 'ocean_cobalt_feisty_forcing_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.thkcello.mat'])

    %% Vertical means
    vSp = mean(nsm,1,'omitnan');
    vSp = squeeze(mean(vSp,2,'omitnan'));
    vSp = squeeze(mean(vSp,2,'omitnan'));

    vLp = mean(nlg,1,'omitnan');
    vLp = squeeze(mean(vLp,2,'omitnan'));
    vLp = squeeze(mean(vLp,2,'omitnan'));

    vSz = mean(nsmz,1,'omitnan');
    vSz = squeeze(mean(vSz,2,'omitnan'));
    vSz = squeeze(mean(vSz,2,'omitnan'));

    vDe = mean(ndet,1,'omitnan');
    vDe = squeeze(mean(vDe,2,'omitnan'));
    vDe = squeeze(mean(vDe,2,'omitnan'));

    vDi = mean(ndi,1,'omitnan');
    vDi = squeeze(mean(vDi,2,'omitnan'));
    vDi = squeeze(mean(vDi,2,'omitnan'));
  
    %% vert sums or means
    iSp = squeeze(sum((nsm.*thkcello),3,'omitnan'));
    iLp = squeeze(sum((nlg.*thkcello),3,'omitnan'));
    iSz = squeeze(sum((nsmz.*thkcello),3,'omitnan'));
    iDe = squeeze(sum((ndet.*thkcello),3,'omitnan'));
    iDi = squeeze(sum((ndi.*thkcello),3,'omitnan'));
    
    %% Time series of vert integral
    tSp = mean(iSp,1,'omitnan');
    tSp = squeeze(mean(tSp,2,'omitnan'));

    tLp = mean(iLp,1,'omitnan');
    tLp = squeeze(mean(tLp,2,'omitnan'));

    tSz = mean(iSz,1,'omitnan');
    tSz = squeeze(mean(tSz,2,'omitnan'));

    tDe = mean(iDe,1,'omitnan');
    tDe = squeeze(mean(tDe,2,'omitnan'));

    tDi = mean(iDi,1,'omitnan');
    tDi = squeeze(mean(tDi,2,'omitnan'));

    %% spatial mean of vert integral
    sSp = mean(iSp,3,'omitnan');
    sLp = mean(iLp,3,'omitnan');
    sSz = mean(iSz,3,'omitnan');
    sDe = mean(iDe,3,'omitnan');
    sDi = mean(iDi,3,'omitnan');
   
    %% put in arrays
    yid = (((y-1)*60)+1):(y*60);
   
    tSP(1,yid) = tSp;
    tLP(1,yid) = tLp;
    tSZ(1,yid) = tSz;
    tDE(1,yid) = tDe;
    tDI(1,yid) = tDi;

    sSP(:,:,y) = sSp;
    sLP(:,:,y) = sLp;
    sSZ(:,:,y) = sSz;
    sDE(:,:,y) = sDe;
    sDI(:,:,y) = sDi;

    vSP(:,y) = vSp;
    vLP(:,y) = vLp;
    vSZ(:,y) = vSz;
    vDE(:,y) = vDe;
    vDI(:,y) = vDi;

end

save([fpath 'ocean_cobalt_tracers_month_z.199001',...
        '-',num2str(en(y)),'12_means.nc'],'tSP','tLP','tSZ','tDE','tDI',...
        'sSP','sLP','sSZ','sDE','sDI',...
        'vSP','vLP','vSZ','vDE','vDI')






