% 0.5 degree global MOM6-COBALTv3 only


clear
close all

%fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';
fpath = '/project/Feisty/NC/Global_COBALT_FEISTY/';

%gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
gpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

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
    load([gpath 'ocean_cobalt_feisty_forcing_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.thkcello.mat'])

    %% Vertical means, area weighted
    [ni,nj] = size(geolat);
    matarea1 = reshape(areacello,ni*nj,1);
    matarea = repmat(matarea1,1,35,60);
    matarea2 = repmat(matarea1,1,60);
    matthk = reshape(thkcello,ni*nj,35,60);

    matSp = reshape(nsm,ni*nj,35,60);
    matLp = reshape(nlg,ni*nj,35,60);
    matSz = reshape(nsmz,ni*nj,35,60);
    matDi = reshape(ndi,ni*nj,35,60);
    matDe = reshape(ndet,ni*nj,35,60);

    %eq
    ear = matarea(eq,:,:);
    eth = matthk(eq,:,:);

    evSp = matSp(eq,:,:);
    evSp = squeeze(sum(evSp.*ear,1,'omitnan') ./ sum(ear,1,'omitnan'));
    evSp = squeeze(mean(evSp,2,'omitnan'));

    evLp = matLp(eq,:,:);
    evLp = squeeze(sum(evLp.*ear,1,'omitnan') ./ sum(ear,1,'omitnan'));
    evLp = squeeze(mean(evLp,2,'omitnan'));

    evSz = matSz(eq,:,:);
    evSz = squeeze(sum(evSz.*ear,1,'omitnan') ./ sum(ear,1,'omitnan'));
    evSz = squeeze(mean(evSz,2,'omitnan'));

    evDi = matDi(eq,:,:);
    evDi = squeeze(sum(evDi.*ear,1,'omitnan') ./ sum(ear,1,'omitnan'));
    evDi = squeeze(mean(evDi,2,'omitnan'));

    evDe = matDe(eq,:,:);
    evDe = squeeze(sum(evDe.*ear,1,'omitnan') ./ sum(ear,1,'omitnan'));
    evDe = squeeze(mean(evDe,2,'omitnan'));


    %subtrop
    sar = matarea(sub,:,:);
    sth = matthk(sub,:,:);

    svSp = matSp(sub,:,:);
    svSp = squeeze(sum(svSp.*sar,1,'omitnan') ./ sum(sar,1,'omitnan'));
    svSp = squeeze(mean(svSp,2,'omitnan'));

    svLp = matLp(sub,:,:);
    svLp = squeeze(sum(svLp.*sar,1,'omitnan') ./ sum(sar,1,'omitnan'));
    svLp = squeeze(mean(svLp,2,'omitnan'));

    svSz = matSz(sub,:,:);
    svSz = squeeze(sum(svSz.*sar,1,'omitnan') ./ sum(sar,1,'omitnan'));
    svSz = squeeze(mean(svSz,2,'omitnan'));

    svDi = matDi(sub,:,:);
    svDi = squeeze(sum(svDi.*sar,1,'omitnan') ./ sum(sar,1,'omitnan'));
    svDi = squeeze(mean(svDi,2,'omitnan'));

    svDe = matDe(sub,:,:);
    svDe = squeeze(sum(svDe.*sar,1,'omitnan') ./ sum(sar,1,'omitnan'));
    svDe = squeeze(mean(svDe,2,'omitnan'));


    %temp/subp
    tar = matarea(tem,:,:);
    tth = matthk(tem,:,:);

    tvSp = matSp(tem,:,:);
    tvSp = squeeze(sum(tvSp.*tar,1,'omitnan') ./ sum(tar,1,'omitnan'));
    tvSp = squeeze(mean(tvSp,2,'omitnan'));

    tvLp = matLp(tem,:,:);
    tvLp = squeeze(sum(tvLp.*tar,1,'omitnan') ./ sum(tar,1,'omitnan'));
    tvLp = squeeze(mean(tvLp,2,'omitnan'));

    tvSz = matSz(tem,:,:);
    tvSz = squeeze(sum(tvSz.*tar,1,'omitnan') ./ sum(tar,1,'omitnan'));
    tvSz = squeeze(mean(tvSz,2,'omitnan'));

    tvDi = matDi(tem,:,:);
    tvDi = squeeze(sum(tvDi.*tar,1,'omitnan') ./ sum(tar,1,'omitnan'));
    tvDi = squeeze(mean(tvDi,2,'omitnan'));

    tvDe = matDe(tem,:,:);
    tvDe = squeeze(sum(tvDe.*tar,1,'omitnan') ./ sum(tar,1,'omitnan'));
    tvDe = squeeze(mean(tvDe,2,'omitnan'));

    %polar
    par = matarea(pol,:,:);
    pth = matthk(pol,:,:);

    pvSp = matSp(pol,:,:);
    pvSp = squeeze(sum(pvSp.*par,1,'omitnan') ./ sum(par,1,'omitnan'));
    pvSp = squeeze(mean(pvSp,2,'omitnan'));

    pvLp = matLp(pol,:,:);
    pvLp = squeeze(sum(pvLp.*par,1,'omitnan') ./ sum(par,1,'omitnan'));
    pvLp = squeeze(mean(pvLp,2,'omitnan'));
    
    pvSz = matSz(pol,:,:);
    pvSz = squeeze(sum(pvSz.*par,1,'omitnan') ./ sum(par,1,'omitnan'));
    pvSz = squeeze(mean(pvSz,2,'omitnan'));

    pvDe = matDe(pol,:,:);
    pvDe = squeeze(sum(pvDe.*par,1,'omitnan') ./ sum(par,1,'omitnan'));
    pvDe = squeeze(mean(pvDe,2,'omitnan'));

    pvDi = matDi(pol,:,:);
    pvDi = squeeze(sum(pvDi.*par,1,'omitnan') ./ sum(par,1,'omitnan'));
    pvDi = squeeze(mean(pvDi,2,'omitnan'));

  
    %% vert & spatial means
    %eq
    emth = mean(eth,3,'omitnan');
    emth = mean(emth,1,'omitnan');
    eiSp = squeeze(sum(evSp.*emth','omitnan'));
    eiLp = squeeze(sum(evLp.*emth','omitnan'));
    eiSz = squeeze(sum(evSz.*emth','omitnan'));
    eiDi = squeeze(sum(evDi.*emth','omitnan'));
    eiDe = squeeze(sum(evDe.*emth','omitnan'));

    %sub
    smth = mean(sth,3,'omitnan');
    smth = mean(smth,1,'omitnan');
    siSp = squeeze(sum(svSp.*smth','omitnan'));
    siLp = squeeze(sum(svLp.*smth','omitnan'));
    siSz = squeeze(sum(svSz.*smth','omitnan'));
    siDi = squeeze(sum(svDi.*smth','omitnan'));
    siDe = squeeze(sum(svDe.*smth','omitnan'));

    %temp
    tmth = mean(tth,3,'omitnan');
    tmth = mean(tmth,1,'omitnan');
    tiSp = squeeze(sum(tvSp.*tmth','omitnan'));
    tiLp = squeeze(sum(tvLp.*tmth','omitnan'));
    tiSz = squeeze(sum(tvSz.*tmth','omitnan'));
    tiDi = squeeze(sum(tvDi.*tmth','omitnan'));
    tiDe = squeeze(sum(tvDe.*tmth','omitnan'));

    %pol
    pmth = mean(pth,3,'omitnan');
    pmth = mean(pmth,1,'omitnan');
    piSp = squeeze(sum(pvSp.*pmth','omitnan'));
    piLp = squeeze(sum(pvLp.*pmth','omitnan'));
    piSz = squeeze(sum(pvSz.*pmth','omitnan'));
    piDi = squeeze(sum(pvDi.*pmth','omitnan'));
    piDe = squeeze(sum(pvDe.*pmth','omitnan'));


    %% Seasonal cycle of vert integral - TO DO
   
    %% put in arrays
    vSP(:,1) = evSp;
    vSP(:,2) = svSp;
    vSP(:,3) = tvSp;
    vSP(:,4) = pvSp;
    
    vLP(:,1) = evLp;
    vLP(:,2) = svLp;
    vLP(:,3) = tvLp;
    vLP(:,4) = pvLp;

    vSZ(:,1) = evSz;
    vSZ(:,2) = svSz;
    vSZ(:,3) = tvSz;
    vSZ(:,4) = pvSz;
    
    vDI(:,1) = evDi;
    vDI(:,2) = svDi;
    vDI(:,3) = tvDi;
    vDI(:,4) = pvDi;

    vDE(:,1) = evDe;
    vDE(:,2) = svDe;
    vDE(:,3) = tvDe;
    vDE(:,4) = pvDe;

    mSP(1,1) = eiSp;
    mSP(2,1) = siSp;
    mSP(3,1) = tiSp;
    mSP(4,1) = piSp;
    
    mLP(1,1) = eiLp;
    mLP(2,1) = siLp;
    mLP(3,1) = tiLp;
    mLP(4,1) = piLp;

    mSZ(1,1) = eiSz;
    mSZ(2,1) = siSz;
    mSZ(3,1) = tiSz;
    mSZ(4,1) = piSz;
    
    mDI(1,1) = eiDi;
    mDI(2,1) = siDi;
    mDI(3,1) = tiDi;
    mDI(4,1) = piDi;
    
    mDE(1,1) = eiDe;
    mDE(2,1) = siDe;
    mDE(3,1) = tiDe;
    mDE(4,1) = piDe;

end

% save([fpath 'ocean_cobalt_tracers_month_z.199001',...
%         '-',num2str(en(y)),'12_means.mat'],'tSP','tLP','tSZ',...
%         'sSP','sLP','sSZ',...
%         'vSP','vLP','vSZ')

save([fpath 'ocean_cobalt_tracers_month_z.199001',...
        '-',num2str(en(y)),'12_means.mat'],...
        'mSP','mLP','mSZ','mDE','mDI',...
        'vSP','vLP','vSZ','vDE','vDI')






