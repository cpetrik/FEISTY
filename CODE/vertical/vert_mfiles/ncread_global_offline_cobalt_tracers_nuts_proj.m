% 0.5 degree global MOM6-COBALTv3 only
% C, N, O vars

clear
close all

%fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
fpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
gpath = fpath;

%%
%load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet','z_l_units','z_l_long_name','z_l')

%%
ncdisp([fpath 'ocean_cobalt_tracers_month_z.199001-199412.no3.nc'])

%%
st = 1990:5:2019;
en = 1994:5:2019;

tNO3 = nan*ones(1,5*length(st));
tNH4 = tNO3;
tCHL = tNO3;
tO2  = tNO3;

vNO3 = nan*ones(35,length(st));
vNH4 = vNO3;
vCHL = vNO3;
vO2  = vNO3;

sNO3 = nan*ones(720,576,length(st));
sNH4 = sNO3;
sCHL = sNO3;
sO2  = sNO3;
sCHLs = sNO3;

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
    
    %% Vertical means
    vNo3 = mean(no3,1,'omitnan');
    vNo3 = squeeze(mean(vNo3,2,'omitnan'));
    vNo3 = squeeze(mean(vNo3,2,'omitnan'));

    vNh4 = mean(nh4,1,'omitnan');
    vNh4 = squeeze(mean(vNh4,2,'omitnan'));
    vNh4 = squeeze(mean(vNh4,2,'omitnan'));

    vChl = mean(chl,1,'omitnan');
    vChl = squeeze(mean(vChl,2,'omitnan'));
    vChl = squeeze(mean(vChl,2,'omitnan'));

    vo2 = mean(o2,1,'omitnan');
    vo2 = squeeze(mean(vo2,2,'omitnan'));
    vo2 = squeeze(mean(vo2,2,'omitnan'));
  
    %% vert sums or means
    mNo3 = squeeze(sum((no3.*thkcello),3,'omitnan')) ./ squeeze(sum(thkcello,3,'omitnan'));
    mNh4 = squeeze(sum((nh4.*thkcello),3,'omitnan')) ./ squeeze(sum(thkcello,3,'omitnan'));
    iChl = squeeze(sum((chl.*thkcello),3,'omitnan'));
    mo2 = squeeze(sum((o2.*thkcello),3,'omitnan')) ./ squeeze(sum(thkcello,3,'omitnan'));
    
    %% Time series of vert integral
    tNo3 = mean(mNo3,1,'omitnan');
    tNo3 = squeeze(mean(tNo3,2,'omitnan'));

    tNh4 = mean(mNh4,1,'omitnan');
    tNh4 = squeeze(mean(tNh4,2,'omitnan'));

    tChl = mean(iChl,1,'omitnan');
    tChl = squeeze(mean(tChl,2,'omitnan'));

    to2 = mean(mo2,1,'omitnan');
    to2 = squeeze(mean(to2,2,'omitnan'));

    %% spatial mean of vert integral
    sNo3 = mean(mNo3,3,'omitnan');
    sNh4 = mean(mNh4,3,'omitnan');
    sChl = mean(iChl,3,'omitnan');
    so2 = mean(mo2,3,'omitnan');
    sChls = mean(schl,3,'omitnan');
   
    %% put in arrays
    yid = (((y-1)*60)+1):(y*60);
   
    tNO3(1,yid) = tNo3;
    tNH4(1,yid) = tNh4;
    tCHL(1,yid) = tChl;
    tO2(1,yid) = to2;

    vNO3(:,y) = vNo3;
    vNH4(:,y) = vNh4;
    vCHL(:,y) = vChl;
    vO2(:,y) = vo2;

    sNO3(:,:,y) = sNo3;
    sNH4(:,:,y) = sNh4;
    sCHL(:,:,y) = sChl;
    sO2(:,:,y) = so2;
    sCHLs(:,:,y) = sChls;

end

save([fpath 'ocean_cobalt_nuts_month_z.199001',...
        '-',num2str(en(y)),'12_means.mat'],'tNO3','tNH4','tCHL','tO2',...
        'sNO3','sNH4','sCHL','sCHLs','sO2',...
        'vNO3','vNH4','vCHL','vO2')






