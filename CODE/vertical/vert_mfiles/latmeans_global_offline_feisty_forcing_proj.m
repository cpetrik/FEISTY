% 0.5 degree global MOM6-COBALT only
% means and seasonal cycle by latitude bands

clear
close all

%%
%fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';
fpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
gpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/COBALTonly/';

%%
load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet','areacello',...
 'z_l_units','z_l_long_name','z_l','geolat')

%%
%ncdisp([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412.nmdz.nc'])

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
% st = 1990:5:2019;
% en = 1994:5:2019;
st = 1990;
en = 1994;

vMZ = nan*ones(35,4);
vLZ = vMZ;

mMZ = nan*ones(4,1);
mLZ = mMZ;
mDE = mMZ;

%%
for y=1%:6

    %% MZ
    ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.nmdz.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 5
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    nmdz(nmdz>1e19) = nan;

    %% LZ
    ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.nlgz.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 5
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    nlgz(nlgz>1e19) = nan;

    %% Det btm
    ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_2d.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.fntot_btm.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 4
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    fntot_btm(fntot_btm>1e19) = nan;

    %% thkcello
    load([gpath 'ocean_cobalt_feisty_forcing_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.thkcello.mat'])

    %% Vertical means, area weighted
    [ni,nj] = size(geolat);
    matarea1 = reshape(areacello,ni*nj,1);
    matarea = repmat(matarea1,1,35,60);
    matarea2 = repmat(matarea1,1,60);
    matthk = reshape(thkcello,ni*nj,35,60);
    matMz = reshape(nmdz,ni*nj,35,60);
    matLz = reshape(nlgz,ni*nj,35,60);
    matDe = reshape(fntot_btm,ni*nj,60);

    %eq
    ear = matarea(eq,:,:);
    eth = matthk(eq,:,:);

    evMz = matMz(eq,:,:);
    evMz = squeeze(sum(evMz.*ear,1,'omitnan') ./ sum(ear,1,'omitnan'));
    evMz = squeeze(mean(evMz,2,'omitnan'));

    evLz = matLz(eq,:,:);
    evLz = squeeze(sum(evLz.*ear,1,'omitnan') ./ sum(ear,1,'omitnan'));
    evLz = squeeze(mean(evLz,2,'omitnan'));


    %subtrop
    sar = matarea(sub,:,:);
    sth = matthk(sub,:,:);

    svMz = matMz(sub,:,:);
    svMz = squeeze(sum(svMz.*sar,1,'omitnan') ./ sum(sar,1,'omitnan'));
    svMz = squeeze(mean(svMz,2,'omitnan'));

    svLz = matLz(sub,:,:);
    svLz = squeeze(sum(svLz.*sar,1,'omitnan') ./ sum(sar,1,'omitnan'));
    svLz = squeeze(mean(svLz,2,'omitnan'));


    %temp/subp
    tar = matarea(tem,:,:);
    tth = matthk(tem,:,:);

    tvMz = matMz(tem,:,:);
    tvMz = squeeze(sum(tvMz.*tar,1,'omitnan') ./ sum(tar,1,'omitnan'));
    tvMz = squeeze(mean(tvMz,2,'omitnan'));

    tvLz = matLz(tem,:,:);
    tvLz = squeeze(sum(tvLz.*tar,1,'omitnan') ./ sum(tar,1,'omitnan'));
    tvLz = squeeze(mean(tvLz,2,'omitnan'));


    %polar
    par = matarea(pol,:,:);
    pth = matthk(pol,:,:);

    pvMz = matMz(pol,:,:);
    pvMz = squeeze(sum(pvMz.*par,1,'omitnan') ./ sum(par,1,'omitnan'));
    pvMz = squeeze(mean(pvMz,2,'omitnan'));

    pvLz = matLz(pol,:,:);
    pvLz = squeeze(sum(pvLz.*par,1,'omitnan') ./ sum(par,1,'omitnan'));
    pvLz = squeeze(mean(pvLz,2,'omitnan'));
    
    %% vert & spatial means
    %eq
    emth = mean(eth,3,'omitnan');
    emth = mean(emth,1,'omitnan');
    eiMZ = squeeze(sum(evMz.*emth','omitnan'));
    eiLZ = squeeze(sum(evLz.*emth','omitnan'));
    eDet = squeeze(sum(matDe(eq,:).*matarea2(eq,:),1,'omitnan') ./ sum(matarea2(eq,:),1,'omitnan'));
    emDet = mean(eDet,'omitnan');

    %sub
    smth = mean(sth,3,'omitnan');
    smth = mean(smth,1,'omitnan');
    siMZ = squeeze(sum(svMz.*smth','omitnan'));
    siLZ = squeeze(sum(svLz.*smth','omitnan'));
    sDet = squeeze(sum(matDe(sub,:).*matarea2(sub,:),1,'omitnan') ./ sum(matarea2(sub,:),1,'omitnan'));
    smDet = mean(sDet,'omitnan');

    %temp
    tmth = mean(tth,3,'omitnan');
    tmth = mean(tmth,1,'omitnan');
    tiMZ = squeeze(sum(tvMz.*tmth','omitnan'));
    tiLZ = squeeze(sum(tvLz.*tmth','omitnan'));
    tDet = squeeze(sum(matDe(tem,:).*matarea2(tem,:),1,'omitnan') ./ sum(matarea2(tem,:),1,'omitnan'));
    tmDet = mean(tDet,'omitnan');

    %pol
    pmth = mean(pth,3,'omitnan');
    pmth = mean(pmth,1,'omitnan');
    piMZ = squeeze(sum(pvMz.*pmth','omitnan'));
    piLZ = squeeze(sum(pvLz.*pmth','omitnan'));
    pDet = squeeze(sum(matDe(pol,:).*matarea2(pol,:),1,'omitnan') ./ sum(matarea2(pol,:),1,'omitnan'));
    pmDet = mean(pDet,'omitnan');

    %% Seasonal cycle of vert integral - TO DO
    %eq
    emar = squeeze(mean(ear,2,'omitnan'));
    
    emMZ = squeeze(sum((matMz(eq,:,:).*matthk(eq,:,:)),2,'omitnan'));
    emMZ = squeeze(sum(emMZ.*emar,1,'omitnan') ./ sum(emar,1,'omitnan'));
    
    emLZ = squeeze(sum((matLz(eq,:,:).*matthk(eq,:,:)),2,'omitnan'));
    emLZ = squeeze(sum(emLZ.*emar,1,'omitnan') ./ sum(emar,1,'omitnan'));

    %subtrop
    smar = squeeze(mean(sar,2,'omitnan'));
    
    smMZ = squeeze(sum((matMz(sub,:,:).*matthk(sub,:,:)),2,'omitnan'));
    smMZ = squeeze(sum(smMZ.*smar,1,'omitnan') ./ sum(smar,1,'omitnan'));
    
    smLZ = squeeze(sum((matLz(sub,:,:).*matthk(sub,:,:)),2,'omitnan'));
    smLZ = squeeze(sum(smLZ.*smar,1,'omitnan') ./ sum(smar,1,'omitnan'));
    
    %temp/subp
    tmar = squeeze(mean(tar,2,'omitnan'));

    tmMZ = squeeze(sum((matMz(tem,:,:).*matthk(tem,:,:)),2,'omitnan'));
    tmMZ = squeeze(sum(tmMZ.*tmar,1,'omitnan') ./ sum(tmar,1,'omitnan'));
   
    tmLZ = squeeze(sum((matLz(tem,:,:).*matthk(tem,:,:)),2,'omitnan'));
    tmLZ = squeeze(sum(tmLZ.*tmar,1,'omitnan') ./ sum(tmar,1,'omitnan'));

    %polar
    pmar = squeeze(mean(par,2,'omitnan'));

    pmMZ = squeeze(sum((matMz(pol,:,:).*matthk(pol,:,:)),2,'omitnan'));
    pmMZ = squeeze(sum(pmMZ.*pmar,1,'omitnan') ./ sum(pmar,1,'omitnan'));
    
    pmLZ = squeeze(sum((matLz(pol,:,:).*matthk(pol,:,:)),2,'omitnan'));
    pmLZ = squeeze(sum(pmLZ.*pmar,1,'omitnan') ./ sum(pmar,1,'omitnan'));

    stMz = nan*ones(12,4);
    stLz = nan*ones(12,4);
    stDet = nan*ones(12,4);
    % NEED TO SHIFT S HEM OR TAKE MEAN SEPARATELY
    for mo = 1:12
        m = mo:12:60;
        
        stMz(mo,1) = mean(emMZ(m),'omitnan');
        stLz(mo,1) = mean(emLZ(m),'omitnan');
        stDet(mo,1) = mean(eDet(m),'omitnan');

        stMz(mo,2) = mean(smMZ(m),'omitnan');
        stLz(mo,2) = mean(smLZ(m),'omitnan');
        stDet(mo,2) = mean(sDet(m),'omitnan');

        stMz(mo,3) = mean(tmMZ(m),'omitnan');
        stLz(mo,3) = mean(tmLZ(m),'omitnan');
        stDet(mo,3) = mean(tDet(m),'omitnan');

        stMz(mo,4) = mean(pmMZ(m),'omitnan');
        stLz(mo,4) = mean(pmLZ(m),'omitnan');
        stDet(mo,4) = mean(pDet(m),'omitnan');
    end

   
    %% put in arrays
    yid = (((y-1)*60)+1):(y*60);
   
    vMZ(:,1) = evMz;
    vMZ(:,2) = svMz;
    vMZ(:,3) = tvMz;
    vMZ(:,4) = pvMz;
    
    vLZ(:,1) = evLz;
    vLZ(:,2) = svLz;
    vLZ(:,3) = tvLz;
    vLZ(:,4) = pvLz;

    mMZ(1,1) = eiMZ;
    mMZ(2,1) = siMZ;
    mMZ(3,1) = tiMZ;
    mMZ(4,1) = piMZ;
    
    mLZ(1,1) = eiLZ;
    mLZ(2,1) = siLZ;
    mLZ(3,1) = tiLZ;
    mLZ(4,1) = piLZ;
    
    mDE(1,1) = emDet;
    mDE(2,1) = smDet;
    mDE(3,1) = tmDet;
    mDE(4,1) = pmDet;

end

save([fpath 'ocean_cobalt_feisty_forcing_z.199001',...
        '-',num2str(en(y)),'12_means_lat.mat'],'stMz','stLz','stDet',...
        'mMZ','mLZ','mDE','vMZ','vLZ')

% save([spath 'ocean_cobalt_feisty_forcing_z.199001',...
%         '-',num2str(en(y)),'12_means_lat.mat'],'stMz','stLz','stDet',...
%         'mMZ','mLZ','mDE','vMZ','vLZ')





