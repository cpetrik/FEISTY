% 0.5 degree global MOM6-COBALTv3-FEISTY online


clear
close all

%fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';
fpath = '/project/Feisty/NC/Global_COBALT_FEISTY/';

%gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
gpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%%
%load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet','z_l_units','z_l_long_name','z_l')

%%
%ncdisp([fpath 'ocean_feisty_tracers_z.199001-199412.Sf_B.nc'])

%%
st = 1990:5:2019;
en = 1994:5:2019;

tSF = nan*ones(1,60*length(st));
tSP = tSF;
tSD = tSF;
tMF = tSF;
tMP = tSF;
tMD = tSF;
tLP = tSF;
tLD = tSF;
tBE = tSF;

sSF = nan*ones(720,576,length(st));
sSP = sSF;
sSD = sSF;
sMF = sSF;
sMP = sSF;
sMD = sSF;
sLP = sSF;
sLD = sSF;
sBE = sSF;

vSF = nan*ones(35,length(st));
vSP = vSF;
vSD = vSF;
vMF = vSF;
vMP = vSF;
vLP = vSF;

%%
for y=1:2%:6

    %% SF
    ncid = netcdf.open([fpath 'ocean_feisty_tracers_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.Sf_B.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 2
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    Sf_B(Sf_B>1e19) = nan;

    %% SP
    ncid = netcdf.open([fpath 'ocean_feisty_tracers_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.Sp_B.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 2
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    Sp_B(Sp_B>1e19) = nan;

    %% SD
    ncid = netcdf.open([fpath 'ocean_feisty_tracers_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.Sd_B.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 2
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    Sd_B(Sd_B>1e19) = nan;

    %% MF
    ncid = netcdf.open([fpath 'ocean_feisty_tracers_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.Mf_B.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 2
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    Mf_B(Mf_B>1e19) = nan;

    %% MP
    ncid = netcdf.open([fpath 'ocean_feisty_tracers_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.Mp_B.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 2
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    Mp_B(Mp_B>1e19) = nan;

    %% MD
    ncid = netcdf.open([fpath 'ocean_feisty_tracers_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.Md_B.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 2
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    Md_B(Md_B>1e19) = nan;

    %% LP
    ncid = netcdf.open([fpath 'ocean_feisty_tracers_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.Lp_B.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 2
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    Lp_B(Lp_B>1e19) = nan;

    % LD
    ncid = netcdf.open([fpath 'ocean_feisty_tracers_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.Ld_B.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 2
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    Ld_B(Ld_B>1e19) = nan;

    % BE
    ncid = netcdf.open([fpath 'ocean_feisty_tracers_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.BE_B.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 2
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    BE_B(BE_B>1e19) = nan;

    %% thkcello
    load([gpath 'ocean_cobalt_feisty_forcing_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.thkcello.mat'])

    %% Vertical means
    vSf = mean(Sf_B,1,'omitnan');
    vSf = squeeze(mean(vSf,2,'omitnan'));
    vSf = squeeze(mean(vSf,2,'omitnan'));

    vSp = mean(Sp_B,1,'omitnan');
    vSp = squeeze(mean(vSp,2,'omitnan'));
    vSp = squeeze(mean(vSp,2,'omitnan'));

    vSd = mean(Sd_B,1,'omitnan');
    vSd = squeeze(mean(vSd,2,'omitnan'));
    vSd = squeeze(mean(vSd,2,'omitnan'));

    vMf = mean(Mf_B,1,'omitnan');
    vMf = squeeze(mean(vMf,2,'omitnan'));
    vMf = squeeze(mean(vMf,2,'omitnan'));

    vMp = mean(Mp_B,1,'omitnan');
    vMp = squeeze(mean(vMp,2,'omitnan'));
    vMp = squeeze(mean(vMp,2,'omitnan'));

    vLp = mean(Lp_B,1,'omitnan');
    vLp = squeeze(mean(vLp,2,'omitnan'));
    vLp = squeeze(mean(vLp,2,'omitnan'));

    vLd = mean(Ld_B,1,'omitnan');
    vLd = squeeze(mean(vLd,2,'omitnan'));
    vLd = squeeze(sum(vLd,2,'omitnan'));
    
  
    %% vert sums
    iSf = squeeze(sum((Sf_B.*thkcello),3,'omitnan'));
    iSp = squeeze(sum((Sp_B.*thkcello),3,'omitnan'));
    iSd = squeeze(sum((Sd_B.*thkcello),3,'omitnan'));
    iMf = squeeze(sum((Mf_B.*thkcello),3,'omitnan'));
    iMp = squeeze(sum((Mp_B.*thkcello),3,'omitnan'));
    iLp = squeeze(sum((Lp_B.*thkcello),3,'omitnan'));
   
    %% Time series (of vert integral)
    tSf = mean(iSf,1,'omitnan');
    tSf = squeeze(mean(tSf,2,'omitnan'));

    tSp = mean(iSp,1,'omitnan');
    tSp = squeeze(mean(tSp,2,'omitnan'));

    tSd = mean(iSd,1,'omitnan');
    tSd = squeeze(mean(tSd,2,'omitnan'));

    tMf = mean(iMf,1,'omitnan');
    tMf = squeeze(mean(tMf,2,'omitnan'));

    tMp = mean(iMp,1,'omitnan');
    tMp = squeeze(mean(tMp,2,'omitnan'));

    tLp = mean(iLp,1,'omitnan');
    tLp = squeeze(mean(tLp,2,'omitnan'));

    tMd = mean(Md_B(:,:,35,:),1,'omitnan');
    tMd = squeeze(mean(tMd,2,'omitnan'));

    tLd = mean(Ld_B(:,:,35,:),1,'omitnan');
    tLd = squeeze(mean(tLd,2,'omitnan'));

    tBe = mean(BE_B(:,:,35,:),1,'omitnan');
    tBe = squeeze(mean(tBe,2,'omitnan'));

    %% spatial mean (of vert integral)
    sSf = mean(iSf,3,'omitnan');
    sSp = mean(iSp,3,'omitnan');
    sSd = mean(iSd,3,'omitnan');
    sMf = mean(iMf,3,'omitnan');
    sMp = mean(iMp,3,'omitnan');
    sLp = mean(iLp,3,'omitnan');

    sMd = mean(Md_B(:,:,35,:),4,'omitnan');
    sLd = mean(Ld_B(:,:,35,:),4,'omitnan');
    sBe = mean(BE_B(:,:,35,:),4,'omitnan');
   
    %% put in arrays
    yid = (((y-1)*60)+1):(y*60);
   
    tSF(1,yid) = tSf;
    tSP(1,yid) = tSp;
    tSD(1,yid) = tSd;
    tMF(1,yid) = tMf;
    tMP(1,yid) = tMp;
    tLP(1,yid) = tLp;
    tMD(1,yid) = tMd';
    tLD(1,yid) = tLd';
    tLP(1,yid) = tBe';

    sSF(:,:,y) = sSf;
    sSP(:,:,y) = sSp;
    sSD(:,:,y) = sSd;
    sMF(:,:,y) = sMf;
    sMP(:,:,y) = sMp;
    sMD(:,:,y) = sMd;
    sLP(:,:,y) = sLp;
    sLD(:,:,y) = sLd;
    sBE(:,:,y) = sBe;

    vSF(:,y) = vSf;
    vSP(:,y) = vSp;
    vSD(:,y) = vSd;
    vMF(:,y) = vMf;
    vMP(:,y) = vMp;
    vLP(:,y) = vLp;

end

save([fpath 'ocean_feisty_tracers_z.199001',...
        '-',num2str(en(y)),'12_means.mat'],'tSF','tSP','tSD',...
        'tMF','tMP','tMD','tLP','tLD','tBE','sSF','sSP','sSD',...
        'sMF','sMP','sMD','sLP','sLD','sBE','vSF','vSP','vSD',...
        'vMF','vMP','vLP')






