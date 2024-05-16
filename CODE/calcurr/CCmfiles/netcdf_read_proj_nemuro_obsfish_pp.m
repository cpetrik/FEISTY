% FEISTY output at all locations

clear
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';
esm1 = {'IPSL','GFDL','HAD'};
esm2 = {'ipsl','gfdl','hadley'};
harv = 'All_fishobs';

%%
for z=1:3

    vers = esm1{z};
    vers2 = esm2{z};

    fpath=['/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/',vers,'/'];
    Cdir = ['/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/',vers,'down/'];

    load([Cdir 'Data_grid_nemuro_',vers2,'.mat'],'GRD')

    %% MF
    ncid = netcdf.open([fpath 'Project_' vers '_' harv '_pp_med_f.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    [nid,nt] = size(popprod);

    popprod(popprod(:)>1) = 0; 
    MF.pp = popprod;
    clear popprod

    % LP
    ncid = netcdf.open([fpath 'Project_' vers '_' harv '_pp_lrg_p.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    popprod(popprod(:)>1) = 0; 
    LP.pp = popprod;
    clear popprod

    % LD
    ncid = netcdf.open([fpath 'Project_' vers '_' harv '_pp_lrg_d.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    popprod(popprod(:)>1) = 0; 
    LD.pp = popprod;
    clear popprod 

    %% Take means
    %Time
    mf_tmpp=mean(MF.pp,1,"omitnan");
    lp_tmpp=mean(LP.pp,1,"omitnan");
    ld_tmpp=mean(LD.pp,1,"omitnan");

    %% All years
    %lyr=time((end-12+1):end);
    %lyr=1:12;
    mf_smpp=mean(MF.pp,2,"omitnan");
    lp_smpp=mean(LP.pp,2,"omitnan");
    ld_smpp=mean(LD.pp,2,"omitnan");

    %% Each year
    a = 1:12:nt; % start of each yr
    b = 12:12:nt; % end of each yr
    ymrB = NaN*ones(length(mf_smpp),(nt/12));
    ymppMF = ymrB;
    ymppLP = ymrB;
    ymppLD = ymrB;
    

    for i = 1:(nt/12)
        ymppMF(:,i) = mean(MF.pp(:,a(i):b(i)),2,"omitnan");
        ymppLP(:,i) = mean(LP.pp(:,a(i):b(i)),2,"omitnan");
        ymppLD(:,i) = mean(LD.pp(:,a(i):b(i)),2,"omitnan");

    end

    %%
    save([fpath 'Means_Project_' vers '_' harv '_pp_' cfile '.mat'],...
        'mf_smpp','lp_smpp','ld_smpp',...
        'mf_tmpp','lp_tmpp','ld_tmpp',...
        'ymppMF','ymppLP','ymppLD');

end

