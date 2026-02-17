% 0.5 degree global MOM6-COBALTv3 only


clear
close all

%fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';
fpath = '/project/Feisty/NC/Global_COBALT_FEISTY/';

%gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
gpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%%
%load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet','z_l_units','z_l_long_name','z_l')

%%
%ncdisp([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412.nmdz.nc'])

%%
st = 1990:5:2019;
en = 1994:5:2019;

tMZ = nan*ones(1,60*length(st));
tMH = tMZ;
tLZ = tMZ;
tLH = tMZ;
tTP = tMZ;
tTB = tMZ;
tDE = tMZ;

sMZ = nan*ones(720,576,length(st));
sMH = sMZ;
sLZ = sMZ;
sLH = sMZ;
sTP = sMZ;
sTB = sMZ;
sDE = sMZ;

vMZ = nan*ones(35,length(st));
vMH = vMZ;
vLZ = vMZ;
vLH = vMZ;

%%
for y=1:2%:6

    %% thetao
    ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.thetao.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 6
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    thetao(thetao>1e19) = nan;

    %% tbtm
    ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_2d.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.tob.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 7
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    tob(tob>1e19) = nan;

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

    %% MZ loss
    ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.jhploss_n_Mdz.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 5
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    jhploss_n_Mdz(jhploss_n_Mdz>1e19) = nan;

    %% LZ loss
    ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_z.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.jhploss_n_Lgz.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 5
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    jhploss_n_Lgz(jhploss_n_Lgz>1e19) = nan;

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

    %% Vertical means
    vMz = mean(nmdz,1,'omitnan');
    vMz = squeeze(mean(vMz,2,'omitnan'));
    vMz = squeeze(mean(vMz,2,'omitnan'));

    vMhp = mean(jhploss_n_Mdz,1,'omitnan');
    vMhp = squeeze(mean(vMhp,2,'omitnan'));
    vMhp = squeeze(mean(vMhp,2,'omitnan'));

    vLz = mean(nlgz,1,'omitnan');
    vLz = squeeze(mean(vLz,2,'omitnan'));
    vLz = squeeze(mean(vLz,2,'omitnan'));

    vLhp = mean(jhploss_n_Lgz,1,'omitnan');
    vLhp = squeeze(mean(vLhp,2,'omitnan'));
    vLhp = squeeze(mean(vLhp,2,'omitnan'));
  
    %% vert sums or means
    iMZ = squeeze(sum((nmdz.*thkcello),3,'omitnan'));
    iMH = squeeze(sum((jhploss_n_Mdz.*thkcello),3,'omitnan'));
    
    iLZ = squeeze(sum((nlgz.*thkcello),3,'omitnan'));
    iLH = squeeze(sum((jhploss_n_Lgz.*thkcello),3,'omitnan'));

    mTP = squeeze(sum((thetao(:,:,1:10,:).*thkcello(:,:,1:10,:)),3,'omitnan') ./ sum(thkcello(:,:,1:10,:),3,'omitnan'));

    %% Time series of vert integral
    tMz = mean(iMZ,1,'omitnan');
    tMz = squeeze(mean(tMz,2,'omitnan'));

    tMhp = mean(iMH,1,'omitnan');
    tMhp = squeeze(mean(tMhp,2,'omitnan'));

    tLz = mean(iLZ,1,'omitnan');
    tLz = squeeze(mean(tLz,2,'omitnan'));

    tLhp = mean(iLH,1,'omitnan');
    tLhp = squeeze(mean(tLhp,2,'omitnan'));

    tTp = mean(mTP,1,'omitnan');
    tTp = squeeze(mean(tTp,2,'omitnan'));

    tTb = mean(tob,1,'omitnan');
    tTb = squeeze(mean(tTb,2,'omitnan'));

    tDet = mean(fntot_btm,1,'omitnan');
    tDet = squeeze(mean(tDet,2,'omitnan'));

    %% spatial mean of vert integral
    sMz = mean(iMZ,3,'omitnan');
    sMhp = mean(iMH,3,'omitnan');
    sLz = mean(iLZ,3,'omitnan');
    sLhp = mean(iLH,3,'omitnan');
    sTp = mean(mTP,3,'omitnan');
    sTb = mean(tob,3,'omitnan');
    sDet = mean(fntot_btm,3,'omitnan');
   
    %% put in arrays
    yid = (((y-1)*60)+1):(y*60);
   
    tMZ(1,yid) = tMz;
    tMH(1,yid) = tMhp;
    tLZ(1,yid) = tLz;
    tLH(1,yid) = tLhp;
    tTP(1,yid) = tTp;
    tTB(1,yid) = tTb';
    tDE(1,yid) = tDet';

    sMZ(:,:,y) = sMz;
    sMH(:,:,y) = sMhp;
    sLZ(:,:,y) = sLz;
    sLH(:,:,y) = sLhp;
    sTP(:,:,y) = sTp;
    sTB(:,:,y) = sTb;
    sDE(:,:,y) = sDet;

    vMZ(:,y) = vMz;
    vMH(:,y) = vMhp;
    vLZ(:,y) = vLz;
    vLH(:,y) = vLhp;

end

save([fpath 'ocean_cobalt_feisty_forcing_z.199001',...
        '-',num2str(en(y)),'12_means.nc'],'tMZ','tMH','tLZ','tLH','tTP',...
        'tTB','tDE','sMZ','sMH','sLZ','sLH','sTP','sTB','sDE',...
        'vMZ','vMH','vLZ','vLH')






