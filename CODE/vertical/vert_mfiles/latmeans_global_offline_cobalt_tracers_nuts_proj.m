% 0.5 degree global MOM6-COBALTv3-FEISTY online
% C, N, O vars

clear
close all

%fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
fpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
gpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%%
load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet','z_l_units','z_l_long_name','z_l')

%%
ncdisp([fpath 'ocean_cobalt_tracers_month_z.199001-199412.no3.nc'])

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

vNO3 = nan*ones(35,4);
vNH4 = vNO3;
vCHL = vNO3;
vO2  = vNO3;

mNO3 = nan*ones(4,1);
mNH4 = mNO3;
mCHL = mNO3;
mO2  = mNO3;
mCHLs = mNO3;

%%
for y=1%:6

    %% NO3 mol/kg
    ncid = netcdf.open([fpath 'ocean_cobalt_tracers_month_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.no3.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 2:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    no3(no3>1e19) = nan;

    %% NH4
    ncid = netcdf.open([fpath 'ocean_cobalt_tracers_month_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.nh4.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 2:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    nh4(nh4>1e19) = nan;

    %% CHL
    ncid = netcdf.open([fpath 'ocean_cobalt_tracers_month_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.chl.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 2:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    chl(chl>1e19) = nan;

    schl = squeeze(chl(:,:,1,:));

    %% O2 'mol m-2' - Annual
    ncid = netcdf.open([fpath 'ocean_cobalt_tracers_month_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.o2.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 2:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    o2(o2>1e19) = nan;

    %% thkcello
    load([gpath 'ocean_cobalt_feisty_forcing_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.thkcello.mat'])
    
    %% Vertical means, area weighted
    [ni,nj] = size(geolat);
    matarea1 = reshape(areacello,ni*nj,1);
    matarea = repmat(matarea1,1,35,60);
    matarea2 = repmat(matarea1,1,60);
    matthk = reshape(thkcello,ni*nj,35,60);
    
    matNo3 = reshape(no3,ni*nj,35,60);
    matNh4 = reshape(nh4,ni*nj,35,60);
    mato2 = reshape(o2,ni*nj,35,60);
    matChl = reshape(chl,ni*nj,35,60);
    matSChl = reshape(schl,ni*nj,60);

    %eq
    ear = matarea(eq,:,:);
    eth = matthk(eq,:,:);

    evNo3 = matNo3(eq,:,:);
    evNo3 = squeeze(sum(evNo3.*ear,1,'omitnan') ./ sum(ear,1,'omitnan'));
    evNo3 = squeeze(mean(evNo3,2,'omitnan'));

    evNh4 = matNh4(eq,:,:);
    evNh4 = squeeze(sum(evNh4.*ear,1,'omitnan') ./ sum(ear,1,'omitnan'));
    evNh4 = squeeze(mean(evNh4,2,'omitnan'));

    evo2 = mato2(eq,:,:);
    evo2 = squeeze(sum(evo2.*ear,1,'omitnan') ./ sum(ear,1,'omitnan'));
    evo2 = squeeze(mean(evo2,2,'omitnan'));

    evChl = matChl(eq,:,:);
    evChl = squeeze(sum(evChl.*ear,1,'omitnan') ./ sum(ear,1,'omitnan'));
    evChl = squeeze(mean(evChl,2,'omitnan'));


    %subtrop
    sar = matarea(sub,:,:);
    sth = matthk(sub,:,:);

    svNo3 = matNo3(sub,:,:);
    svNo3 = squeeze(sum(svNo3.*sar,1,'omitnan') ./ sum(sar,1,'omitnan'));
    svNo3 = squeeze(mean(svNo3,2,'omitnan'));

    svNh4 = matNh4(sub,:,:);
    svNh4 = squeeze(sum(svNh4.*sar,1,'omitnan') ./ sum(sar,1,'omitnan'));
    svNh4 = squeeze(mean(svNh4,2,'omitnan'));

    svo2 = mato2(sub,:,:);
    svo2 = squeeze(sum(svo2.*sar,1,'omitnan') ./ sum(sar,1,'omitnan'));
    svo2 = squeeze(mean(svo2,2,'omitnan'));

    svChl = matChl(sub,:,:);
    svChl = squeeze(sum(svChl.*sar,1,'omitnan') ./ sum(sar,1,'omitnan'));
    svChl = squeeze(mean(svChl,2,'omitnan'));


    %temp/subp
    tar = matarea(tem,:,:);
    tth = matthk(tem,:,:);

    tvNo3 = matNo3(tem,:,:);
    tvNo3 = squeeze(sum(tvNo3.*tar,1,'omitnan') ./ sum(tar,1,'omitnan'));
    tvNo3 = squeeze(mean(tvNo3,2,'omitnan'));

    tvNh4 = matNh4(tem,:,:);
    tvNh4 = squeeze(sum(tvNh4.*tar,1,'omitnan') ./ sum(tar,1,'omitnan'));
    tvNh4 = squeeze(mean(tvNh4,2,'omitnan'));

    tvo2 = mato2(tem,:,:);
    tvo2 = squeeze(sum(tvo2.*tar,1,'omitnan') ./ sum(tar,1,'omitnan'));
    tvo2 = squeeze(mean(tvo2,2,'omitnan'));

    tvChl = matChl(tem,:,:);
    tvChl = squeeze(sum(tvChl.*tar,1,'omitnan') ./ sum(tar,1,'omitnan'));
    tvChl = squeeze(mean(tvChl,2,'omitnan'));


    %polar
    par = matarea(pol,:,:);
    pth = matthk(pol,:,:);

    pvNo3 = matNo3(pol,:,:);
    pvNo3 = squeeze(sum(pvNo3.*par,1,'omitnan') ./ sum(par,1,'omitnan'));
    pvNo3 = squeeze(mean(pvNo3,2,'omitnan'));

    pvNh4 = matNh4(pol,:,:);
    pvNh4 = squeeze(sum(pvNh4.*par,1,'omitnan') ./ sum(par,1,'omitnan'));
    pvNh4 = squeeze(mean(pvNh4,2,'omitnan'));
    
    pvo2 = mato2(pol,:,:);
    pvo2 = squeeze(sum(pvo2.*par,1,'omitnan') ./ sum(par,1,'omitnan'));
    pvo2 = squeeze(mean(pvo2,2,'omitnan'));

    pvChl = matChl(pol,:,:);
    pvChl = squeeze(sum(pvChl.*par,1,'omitnan') ./ sum(par,1,'omitnan'));
    pvChl = squeeze(mean(pvChl,2,'omitnan'));

  
    %% vert & spatial means
    %eq
    emth = mean(eth,3,'omitnan');
    emth = mean(emth,1,'omitnan');
    eiNo3 = squeeze(sum(evNo3.*emth','omitnan'));
    eiNh4 = squeeze(sum(evNh4.*emth','omitnan'));
    eio2 = squeeze(sum(evo2.*emth','omitnan'));
    eiChl = squeeze(sum(evChl.*emth','omitnan'));
    eSChl = squeeze(sum(matSChl(eq,:).*matarea2(eq,:),1,'omitnan') ./ sum(matarea2(eq,:),1,'omitnan'));
    emSChl = mean(eSChl,'omitnan');

    %sub
    smth = mean(sth,3,'omitnan');
    smth = mean(smth,1,'omitnan');
    siNo3 = squeeze(sum(svNo3.*smth','omitnan'));
    siNh4 = squeeze(sum(svNh4.*smth','omitnan'));
    sio2 = squeeze(sum(svo2.*smth','omitnan'));
    siChl = squeeze(sum(svChl.*smth','omitnan'));
    sSChl = squeeze(sum(matSChl(sub,:).*matarea2(sub,:),1,'omitnan') ./ sum(matarea2(sub,:),1,'omitnan'));
    smSChl = mean(sSChl,'omitnan');

    %temp
    tmth = mean(tth,3,'omitnan');
    tmth = mean(tmth,1,'omitnan');
    tiNo3 = squeeze(sum(tvNo3.*tmth','omitnan'));
    tiNh4 = squeeze(sum(tvNh4.*tmth','omitnan'));
    tio2 = squeeze(sum(tvo2.*tmth','omitnan'));
    tiChl = squeeze(sum(tvChl.*tmth','omitnan'));
    tSChl = squeeze(sum(matSChl(tem,:).*matarea2(tem,:),1,'omitnan') ./ sum(matarea2(tem,:),1,'omitnan'));
    tmSChl = mean(tSChl,'omitnan');

    %pol
    pmth = mean(pth,3,'omitnan');
    pmth = mean(pmth,1,'omitnan');
    piNo3 = squeeze(sum(pvNo3.*pmth','omitnan'));
    piNh4 = squeeze(sum(pvNh4.*pmth','omitnan'));
    pio2 = squeeze(sum(pvo2.*pmth','omitnan'));
    piChl = squeeze(sum(pvChl.*pmth','omitnan'));
    pSChl = squeeze(sum(matSChl(pol,:).*matarea2(pol,:),1,'omitnan') ./ sum(matarea2(pol,:),1,'omitnan'));
    pmSChl = mean(pSChl,'omitnan');


    %% Seasonal cycle of vert integral - TO DO
   
    %% put in arrays
    vNO3(:,1) = evNo3;
    vNO3(:,2) = svNo3;
    vNO3(:,3) = tvNo3;
    vNO3(:,4) = pvNo3;
    
    vNH4(:,1) = evNh4;
    vNH4(:,2) = svNh4;
    vNH4(:,3) = tvNh4;
    vNH4(:,4) = pvNh4;

    vO2(:,1) = evo2;
    vO2(:,2) = svo2;
    vO2(:,3) = tvo2;
    vO2(:,4) = pvo2;
    
    vCHL(:,1) = evChl;
    vCHL(:,2) = svChl;
    vCHL(:,3) = tvChl;
    vCHL(:,4) = pvChl;

    mNO3(1,1) = eiNo3;
    mNO3(2,1) = siNo3;
    mNO3(3,1) = tiNo3;
    mNO3(4,1) = piNo3;
    
    mNH4(1,1) = eiNh4;
    mNH4(2,1) = siNh4;
    mNH4(3,1) = tiNh4;
    mNH4(4,1) = piNh4;

    mO2(1,1) = eio2;
    mO2(2,1) = sio2;
    mO2(3,1) = tio2;
    mO2(4,1) = pio2;
    
    mCHL(1,1) = eiChl;
    mCHL(2,1) = siChl;
    mCHL(3,1) = tiChl;
    mCHL(4,1) = piChl;
    
    mSCHL(1,1) = emSChl;
    mSCHL(2,1) = smSChl;
    mSCHL(3,1) = tmSChl;
    mCHLs(4,1) = pmSChl;


end

save([fpath 'ocean_cobalt_nuts_month_z.199001',...
        '-',num2str(en(y)),'12_means_lat.mat'],'tNO3','tNH4','tCHL','tO2',...
        'mNO3','mNH4','mCHL','mCHLs','mO2',...
        'vNO3','vNH4','vCHL','vO2')






