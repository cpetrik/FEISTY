% 0.5 degree global MOM6-COBALTv3-FEISTY only
% C, N, O vars

clear
close all

%fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';
fpath = '/project/Feisty/NC/Global_COBALT_FEISTY/';

%gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
gpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%%
%load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet','z_l_units','z_l_long_name','z_l')

%%
ncdisp([fpath 'ocean_cobalt_tracers_instant.199001-199412.wc_vert_int_poc.nc'])

%%
st = 1990:5:2019;
en = 1994:5:2019;

tPOC = nan*ones(1,5*length(st));
tDOC = tPOC;
tDIC = tPOC;
tO2  = tPOC;
tFN = nan*ones(1,60*length(st));

sPOC = nan*ones(720,576,length(st));
sDOC = sPOC;
sDIC = sPOC;
sO2  = sPOC;
sFN  = sPOC;

%%
for y=1:2%:6

    %% Total particulate organic carbon vertical integral 'mol m-2' - Annual
    ncid = netcdf.open([fpath 'ocean_cobalt_tracers_instant.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.wc_vert_int_poc.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    wc_vert_int_poc(wc_vert_int_poc>1e19) = nan;

    %% DOC  - Annual
    ncid = netcdf.open([fpath 'ocean_cobalt_tracers_instant.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.wc_vert_int_doc.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    wc_vert_int_doc(wc_vert_int_doc>1e19) = nan;

    %% DIC  - Annual
    ncid = netcdf.open([fpath 'ocean_cobalt_tracers_instant.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.wc_vert_int_dic.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    wc_vert_int_dic(wc_vert_int_dic>1e19) = nan;

    %% O2 'mol m-2' - Annual
    ncid = netcdf.open([fpath 'ocean_cobalt_tracers_instant.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.wc_vert_int_o2.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    wc_vert_int_o2(wc_vert_int_o2>1e19) = nan;

    %% Total N sinking flux 'mol m-2 s-1' - Monthly
    ncid = netcdf.open([fpath 'ocean_cobalt_fdet_100.',...
        num2str(st(y)),'01-',num2str(en(y)),'12.fntot_100.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1e20) = NaN;']);
    end
    netcdf.close(ncid);

    fntot_100(fntot_100>1e19) = nan;
    
    %% Time series of vert integral
    tPoc = mean(wc_vert_int_poc,1,'omitnan');
    tPoc = squeeze(mean(tPoc,2,'omitnan'));

    tDoc = mean(wc_vert_int_doc,1,'omitnan');
    tDoc = squeeze(mean(tDoc,2,'omitnan'));

    tDic = mean(wc_vert_int_dic,1,'omitnan');
    tDic = squeeze(mean(tDic,2,'omitnan'));

    to2 = mean(wc_vert_int_o2,1,'omitnan');
    to2 = squeeze(mean(to2,2,'omitnan'));

    tFn = mean(fntot_100,1,'omitnan');
    tFn = squeeze(mean(tFn,2,'omitnan'));

    %% spatial mean of vert integral
    sPoc = mean(wc_vert_int_poc,3,'omitnan');
    sDoc = mean(wc_vert_int_doc,3,'omitnan');
    sDic = mean(wc_vert_int_dic,3,'omitnan');
    so2 = mean(wc_vert_int_o2,3,'omitnan');
    sFn = mean(fntot_100,3,'omitnan');
   
    %% put in arrays
    yid = (((y-1)*5)+1):(y*5);
    mid = (((y-1)*60)+1):(y*60);
   
    tPOC(1,yid) = tPoc;
    tDOC(1,yid) = tDoc;
    tDIC(1,yid) = tDic;
    tO2(1,yid) = to2;
    tFN(1,mid) = tFn;

    sPOC(:,:,y) = sPoc;
    sDOC(:,:,y) = sDoc;
    sDIC(:,:,y) = sDic;
    sO2(:,:,y) = so2;
    sFN(:,:,y) = sFn;

end

save([fpath 'ocean_cobalt_tracers_instant.199001',...
        '-',num2str(en(y)),'12_means.nc'],'tPOC','tFN','tDOC','tDIC','tO2',...
        'sPOC','sFN','sDOC','sDIC','sO2')






