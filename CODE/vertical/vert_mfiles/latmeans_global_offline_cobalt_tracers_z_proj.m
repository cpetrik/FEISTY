% 0.5 degree global MOM6-COBALTv3 only


clear
close all

%fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
fpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%%
load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet','areacello',...
 'z_l_units','z_l_long_name','z_l','geolat')

%%
%ncdisp([fpath 'ocean_cobalt_tracers_month_z.199001-199412.nlg.nc'])

%% lat bands
eq = find(geolat(:)<=10 & geolat(:)>=-10);
elat = geolat(eq);

sub = find(abs(geolat(:))<=40 & abs(geolat(:))>10);
slat = geolat(sub);

tem = find(abs(geolat(:))<=60 & abs(geolat(:))>40);
tlat = geolat(tem);

pol = find(geolat(:)>60 | geolat(:)<-60);
plat = geolat(pol);

figure
subplot(2,2,1)
histogram(elat(:))
subplot(2,2,2)
histogram(slat(:))
subplot(2,2,3)
histogram(tlat(:))
subplot(2,2,4)
histogram(plat(:))

%%
st = 1990;%:5:2019;
en = 1994;%:5:2019;

mSP = nan*ones(4,1);
mLP = mSP;
mSZ = mSP;
mDE = mSP;
mDI = mSP;

vSP = nan*ones(35,4);
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
    load([fpath 'ocean_cobalt_feisty_forcing_z.',...
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

    mSP(:,:,y) = sSp;
    mLP(:,:,y) = sLp;
    mSZ(:,:,y) = sSz;
    mDE(:,:,y) = sDe;
    mDI(:,:,y) = sDi;

    vSP(:,y) = vSp;
    vLP(:,y) = vLp;
    vSZ(:,y) = vSz;
    vDE(:,y) = vDe;
    vDI(:,y) = vDi;

end

% save([fpath 'ocean_cobalt_tracers_month_z.199001',...
%         '-',num2str(en(y)),'12_means.mat'],'tSP','tLP','tSZ',...
%         'sSP','sLP','sSZ',...
%         'vSP','vLP','vSZ')

save([fpath 'ocean_cobalt_tracers_month_z.199001',...
        '-',num2str(en(y)),'12_means.mat'],'tSP','tLP','tSZ','tDE','tDI',...
        'mSP','mLP','mSZ','mDE','mDI',...
        'vSP','vLP','vSZ','vDE','vDI')






