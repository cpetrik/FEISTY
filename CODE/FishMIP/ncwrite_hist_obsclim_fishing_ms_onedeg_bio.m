% FEISTY output at all locations
% Years 1960-2010 will be good
% annual means for biomass and annual sums for catches are okay

clear
close all

%%
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];

mods = {'All_fishobs_assessment','All_fishobs_effective','All_fishobs_nominal',...
    'All_fishobs_FFmsy_creep','All_fishobs_FFmsy_nominal','All_fishobs_FFmsymax_creep',...
    'All_fishobs_FFmsymax_nominal','All_fishobs_FFmsymin_creep','All_fishobs_FFmsymin_nominal'};

%% Grid info
cpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';
load([cpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat']);
load([cpath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat'], 'GRD');
load([cpath 'cellarea_onedeg.mat']);

ID = GRD.ID;

%%
for i=4:length(mods)
    alt = mods{i};

    %% annual mean biomasses
    load([fpath 'Means_Hist_obsclim_',alt,'_',cfile,'.mat'],'time',...
        'sf_mean','sp_mean','sd_mean',...
        'mf_mean','mp_mean','md_mean',...
        'lp_mean','ld_mean','b_mean');

    %% annual total catches
    load([fpath 'Catch_anntots_Hist_obsclim_',alt,'_',cfile,'.mat'],...
        'mf_tac','mp_tac','md_tac','lp_tac','ld_tac');

    %% Netcdf OUTPUTS =================================================

    % Map data
    [ni,nj] = size(LON);
    yr = 1961:2010;
    nt = length(yr);

    allFB = sf_mean + mf_mean;
    allPB = sp_mean + mp_mean + lp_mean;
    allDB = sd_mean + md_mean + ld_mean;
    allFC = mf_tac;
    allPC = mp_tac + lp_tac;
    allDC = md_tac + ld_tac;

    %% Reshape to lat,lon,yr
    AllFB = NaN*ones(ni,nj,nt);
    AllPB = NaN*ones(ni,nj,nt);
    AllDB = NaN*ones(ni,nj,nt);
    AllBB = NaN*ones(ni,nj,nt);
    AllFC = NaN*ones(ni,nj,nt);
    AllPC = NaN*ones(ni,nj,nt);
    AllDC = NaN*ones(ni,nj,nt);

    for z=1:nt
        Zf=NaN*ones(ni,nj);
        Zp=NaN*ones(ni,nj);
        Zd=NaN*ones(ni,nj);
        Zb=NaN*ones(ni,nj);
        Cf=NaN*ones(ni,nj);
        Cp=NaN*ones(ni,nj);
        Cd=NaN*ones(ni,nj);

        Zf(GRD.ID)=allFB(:,z);
        Zp(GRD.ID)=allPB(:,z);
        Zd(GRD.ID)=allDB(:,z);
        Zb(GRD.ID)=b_mean(:,z);
        Cf(GRD.ID)=allFC(:,z);
        Cp(GRD.ID)=allPC(:,z);
        Cd(GRD.ID)=allDC(:,z);

        AllFB(:,:,z) = Zf;
        AllPB(:,:,z) = Zp;
        AllDB(:,:,z) = Zd;
        AllBB(:,:,z) = Zb;
        AllFC(:,:,z) = Cf;
        AllPC(:,:,z) = Cp;
        AllDC(:,:,z) = Cd;
    end

    %% negs to zero
    AllFB(AllFB(:)<0) = 0.0;
    AllPB(AllPB(:)<0) = 0.0;
    AllDB(AllDB(:)<0) = 0.0;
    AllBB(AllBB(:)<0) = 0.0;
    AllFC(AllFC(:)<0) = 0.0;
    AllPC(AllPC(:)<0) = 0.0;
    AllDC(AllDC(:)<0) = 0.0;

    %% Quick look
    pb = AllPC(:,:,50);
    db = AllDC(:,:,50);
    cb = AllFC(:,:,50);

    figure(1)
    pcolor(log10(pb + eps))
    shading flat
    colormap('jet')
    colorbar
    caxis([-1 10])
    title('allP')

    figure(2)
    pcolor(log10(db + eps))
    shading flat
    colormap('jet')
    colorbar
    caxis([-1 10])
    title('allD')

    figure(3)
    pcolor(log10(cb + eps))
    shading flat
    colormap('jet')
    colorbar
    caxis([-1 10])
    title('all F')

    %% netcdf write
    % nans to a large number
    AllFB(isnan(AllFB)) = 1.000000020040877e20;
    AllPB(isnan(AllPB)) = 1.000000020040877e20;
    AllDB(isnan(AllDB)) = 1.000000020040877e20;
    AllBB(isnan(AllBB)) = 1.000000020040877e20;
    AllFC(isnan(AllFC)) = 1.000000020040877e20;
    AllPC(isnan(AllPC)) = 1.000000020040877e20;
    AllDC(isnan(AllDC)) = 1.000000020040877e20;

    %%
    close all

    %% Setup netcdf path to store to
    fname1 = ['FEISTY_obsclim_',alt,'_onedeg_1961-2010_'];
    fname3 = '.nc';

    file_tfb = [fpath fname1 'tfb' fname3];
    file_tpb = [fpath fname1 'tpb' fname3];
    file_tdb = [fpath fname1 'tdb' fname3];
    file_tbb = [fpath fname1 'tbb' fname3];
    file_tfc = [fpath fname1 'tfc' fname3];
    file_tpc = [fpath fname1 'tpc' fname3];
    file_tdc = [fpath fname1 'tdc' fname3];

    [ni,nj,nt] = size(AllPC);

    %Use Netcdf4 classic
    cmode = netcdf.getConstant('NETCDF4');
    cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

    %% tfb
    ncidFB = netcdf.create(file_tfb,cmode);

    time_dim = netcdf.defDim(ncidFB,'time',nt);
    lon_dim = netcdf.defDim(ncidFB,'nlon',ni);
    lat_dim = netcdf.defDim(ncidFB,'nlat',nj);

    vidtFB = netcdf.defVar(ncidFB,'time','NC_DOUBLE',time_dim);
    netcdf.putAtt(ncidFB,vidtFB,'long_name','time');
    netcdf.putAtt(ncidFB,vidtFB,'standard_name','time');
    netcdf.putAtt(ncidFB,vidtFB,'calendar','365_day');
    netcdf.putAtt(ncidFB,vidtFB,'axis','T');
    netcdf.putAtt(ncidFB,vidtFB,'units','year' );

    vidlon = netcdf.defVar(ncidFB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
    netcdf.putAtt(ncidFB,vidlon,'long_name','longitude');
    netcdf.putAtt(ncidFB,vidlon,'standard_name','longitude');
    netcdf.putAtt(ncidFB,vidlon,'axis','X');

    vidlat = netcdf.defVar(ncidFB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
    netcdf.putAtt(ncidFB,vidlat,'long_name','latitude');
    netcdf.putAtt(ncidFB,vidlat,'standard_name','latitude');
    netcdf.putAtt(ncidFB,vidlat,'axis','Y');

    vidbioFB = netcdf.defVar(ncidFB,'tfc','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
    netcdf.defVarChunking(ncidFB,vidbioFB,'CHUNKED',[10, 10, 1]);
    netcdf.putAtt(ncidFB,vidbioFB,'long_name','Mean biomass of Forage Fish');
    netcdf.putAtt(ncidFB,vidbioFB,'units','gWW m-2' );
    netcdf.defVarFill(ncidFB,vidbioFB,false,1.000000020040877e20);
    netcdf.putAtt(ncidFB,vidbioFB,'missing value',1.000000020040877e20);

    varid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncidFB,varid,'creation_date',datestr(now));
    netcdf.putAtt(ncidFB,varid,'_FillValue',1.000000020040877e20);
    netcdf.putAtt(ncidFB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
    netcdf.putAtt(ncidFB,varid,'institution','UCSD');
    netcdf.putAtt(ncidFB,varid,'wet weight:C ratio','9:1');

    netcdf.endDef(ncidFB);

    %tfb = single(tfb);
    netcdf.putVar(ncidFB,vidlat,LAT);
    netcdf.putVar(ncidFB,vidlon,LON);
    netcdf.putVar(ncidFB,vidbioFB,AllFB);
    netcdf.putVar(ncidFB,vidtFB,yr);

    netcdf.close(ncidFB);

    %%
    ncdisp(file_tfb)

    %% tpb
    ncidPB = netcdf.create(file_tpb,cmode);

    time_dim = netcdf.defDim(ncidPB,'time',nt);
    lon_dim = netcdf.defDim(ncidPB,'nlon',ni);
    lat_dim = netcdf.defDim(ncidPB,'nlat',nj);

    vidtPB = netcdf.defVar(ncidPB,'time','NC_DOUBLE',time_dim);
    netcdf.putAtt(ncidPB,vidtPB,'long_name','time');
    netcdf.putAtt(ncidPB,vidtPB,'standard_name','time');
    netcdf.putAtt(ncidPB,vidtPB,'units','year' );
    netcdf.putAtt(ncidPB,vidtPB,'calendar','365_day');
    netcdf.putAtt(ncidPB,vidtPB,'axis','T');

    vidlon = netcdf.defVar(ncidPB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
    netcdf.putAtt(ncidPB,vidlon,'long_name','longitude');
    netcdf.putAtt(ncidPB,vidlon,'standard_name','longitude');
    netcdf.putAtt(ncidPB,vidlon,'axis','X');

    vidlat = netcdf.defVar(ncidPB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
    netcdf.putAtt(ncidPB,vidlat,'long_name','latitude');
    netcdf.putAtt(ncidPB,vidlat,'standard_name','latitude');
    netcdf.putAtt(ncidPB,vidlat,'axis','Y');

    vidbioPB = netcdf.defVar(ncidPB,'tpc','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
    netcdf.defVarChunking(ncidPB,vidbioPB,'CHUNKED',[10, 10, 1]);
    netcdf.putAtt(ncidPB,vidbioPB,'long_name','Mean biomass of Large Pelagic Fish');
    netcdf.putAtt(ncidPB,vidbioPB,'units','gWW m-2' );
    netcdf.defVarFill(ncidPB,vidbioPB,false,1.000000020040877e20);
    netcdf.putAtt(ncidPB,vidbioPB,'missing value',1.000000020040877e20);

    varid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncidPB,varid,'creation_date',datestr(now));
    netcdf.putAtt(ncidPB,varid,'_FillValue',1.000000020040877e20);
    netcdf.putAtt(ncidPB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
    netcdf.putAtt(ncidPB,varid,'institution','UCSD');
    netcdf.putAtt(ncidPB,varid,'wet weight:C ratio','9:1');

    netcdf.endDef(ncidPB);

    %tpb = single(tpb);
    netcdf.putVar(ncidPB,vidlat,LAT);
    netcdf.putVar(ncidPB,vidlon,LON);
    netcdf.putVar(ncidPB,vidbioPB,AllPB);
    netcdf.putVar(ncidPB,vidtPB,yr);

    netcdf.close(ncidPB);

    %% tdb
    ncidDB = netcdf.create(file_tdb,cmode);

    time_dim = netcdf.defDim(ncidDB,'time',nt);
    lon_dim = netcdf.defDim(ncidDB,'nlon',ni);
    lat_dim = netcdf.defDim(ncidDB,'nlat',nj);

    vidtDB = netcdf.defVar(ncidDB,'time','NC_DOUBLE',time_dim);
    netcdf.putAtt(ncidDB,vidtDB,'long_name','time');
    netcdf.putAtt(ncidDB,vidtDB,'standard_name','time');
    netcdf.putAtt(ncidDB,vidtDB,'calendar','365_day');
    netcdf.putAtt(ncidDB,vidtDB,'axis','T');
    netcdf.putAtt(ncidDB,vidtDB,'units','year' );

    vidlon = netcdf.defVar(ncidDB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
    netcdf.putAtt(ncidDB,vidlon,'long_name','longitude');
    netcdf.putAtt(ncidDB,vidlon,'standard_name','longitude');
    netcdf.putAtt(ncidDB,vidlon,'axis','X');

    vidlat = netcdf.defVar(ncidDB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
    netcdf.putAtt(ncidDB,vidlat,'long_name','latitude');
    netcdf.putAtt(ncidDB,vidlat,'standard_name','latitude');
    netcdf.putAtt(ncidDB,vidlat,'axis','Y');

    vidbioDB = netcdf.defVar(ncidDB,'tdc','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
    netcdf.defVarChunking(ncidDB,vidbioDB,'CHUNKED',[10, 10, 1]);
    netcdf.putAtt(ncidDB,vidbioDB,'long_name','Mean biomass of Demersal Fish');
    netcdf.putAtt(ncidDB,vidbioDB,'units','gWW m-2' );
    netcdf.defVarFill(ncidDB,vidbioDB,false,1.000000020040877e20);
    netcdf.putAtt(ncidDB,vidbioDB,'missing value',1.000000020040877e20);

    varid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncidDB,varid,'creation_date',datestr(now));
    netcdf.putAtt(ncidDB,varid,'_FillValue',1.000000020040877e20);
    netcdf.putAtt(ncidDB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
    netcdf.putAtt(ncidDB,varid,'institution','UCSD');
    netcdf.putAtt(ncidDB,varid,'wet weight:C ratio','9:1');

    netcdf.endDef(ncidDB);

    %tdb = single(tdb);
    netcdf.putVar(ncidDB,vidlat,LAT);
    netcdf.putVar(ncidDB,vidlon,LON);
    netcdf.putVar(ncidDB,vidbioDB,AllDB);
    netcdf.putVar(ncidDB,vidtDB,yr);

    netcdf.close(ncidDB);

    %% tbb
    ncidBB = netcdf.create(file_tbb,cmode);

    time_dim = netcdf.defDim(ncidBB,'time',nt);
    lon_dim = netcdf.defDim(ncidBB,'nlon',ni);
    lat_dim = netcdf.defDim(ncidBB,'nlat',nj);

    vidtBB = netcdf.defVar(ncidBB,'time','NC_DOUBLE',time_dim);
    netcdf.putAtt(ncidBB,vidtBB,'long_name','time');
    netcdf.putAtt(ncidBB,vidtBB,'standard_name','time');
    netcdf.putAtt(ncidBB,vidtBB,'calendar','365_day');
    netcdf.putAtt(ncidBB,vidtBB,'axis','T');
    netcdf.putAtt(ncidBB,vidtBB,'units','year' );

    vidlon = netcdf.defVar(ncidBB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
    netcdf.putAtt(ncidBB,vidlon,'long_name','longitude');
    netcdf.putAtt(ncidBB,vidlon,'standard_name','longitude');
    netcdf.putAtt(ncidBB,vidlon,'axis','X');

    vidlat = netcdf.defVar(ncidBB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
    netcdf.putAtt(ncidBB,vidlat,'long_name','latitude');
    netcdf.putAtt(ncidBB,vidlat,'standard_name','latitude');
    netcdf.putAtt(ncidBB,vidlat,'axis','Y');

    vidbioBB = netcdf.defVar(ncidBB,'tdc','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
    netcdf.defVarChunking(ncidBB,vidbioBB,'CHUNKED',[10, 10, 1]);
    netcdf.putAtt(ncidBB,vidbioBB,'long_name','Mean biomass of Benthic Inverts');
    netcdf.putAtt(ncidBB,vidbioBB,'units','gWW m-2' );
    netcdf.defVarFill(ncidBB,vidbioBB,false,1.000000020040877e20);
    netcdf.putAtt(ncidBB,vidbioBB,'missing value',1.000000020040877e20);

    varid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncidBB,varid,'creation_date',datestr(now));
    netcdf.putAtt(ncidBB,varid,'_FillValue',1.000000020040877e20);
    netcdf.putAtt(ncidBB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
    netcdf.putAtt(ncidBB,varid,'institution','UCSD');
    netcdf.putAtt(ncidBB,varid,'wet weight:C ratio','9:1');

    netcdf.endDef(ncidBB);

    netcdf.putVar(ncidBB,vidlat,LAT);
    netcdf.putVar(ncidBB,vidlon,LON);
    netcdf.putVar(ncidBB,vidbioBB,AllBB);
    netcdf.putVar(ncidBB,vidtBB,yr);

    netcdf.close(ncidBB);

    %% tfc
    ncidFC = netcdf.create(file_tfc,cmode);

    time_dim = netcdf.defDim(ncidFC,'time',nt);
    lon_dim = netcdf.defDim(ncidFC,'nlon',ni);
    lat_dim = netcdf.defDim(ncidFC,'nlat',nj);

    vidtFC = netcdf.defVar(ncidFC,'time','NC_DOUBLE',time_dim);
    netcdf.putAtt(ncidFC,vidtFC,'long_name','time');
    netcdf.putAtt(ncidFC,vidtFC,'standard_name','time');
    netcdf.putAtt(ncidFC,vidtFC,'calendar','365_day');
    netcdf.putAtt(ncidFC,vidtFC,'axis','T');
    netcdf.putAtt(ncidFC,vidtFC,'units','year' );

    vidlon = netcdf.defVar(ncidFC,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
    netcdf.putAtt(ncidFC,vidlon,'long_name','longitude');
    netcdf.putAtt(ncidFC,vidlon,'standard_name','longitude');
    netcdf.putAtt(ncidFC,vidlon,'axis','X');

    vidlat = netcdf.defVar(ncidFC,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
    netcdf.putAtt(ncidFC,vidlat,'long_name','latitude');
    netcdf.putAtt(ncidFC,vidlat,'standard_name','latitude');
    netcdf.putAtt(ncidFC,vidlat,'axis','Y');

    vidbioFC = netcdf.defVar(ncidFC,'tfc','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
    netcdf.defVarChunking(ncidFC,vidbioFC,'CHUNKED',[10, 10, 1]);
    netcdf.putAtt(ncidFC,vidbioFC,'long_name','Total Annual Catch of Forage Fish');
    netcdf.putAtt(ncidFC,vidbioFC,'units','gWW' );
    netcdf.defVarFill(ncidFC,vidbioFC,false,1.000000020040877e20);
    netcdf.putAtt(ncidFC,vidbioFC,'missing value',1.000000020040877e20);

    varid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncidFC,varid,'creation_date',datestr(now));
    netcdf.putAtt(ncidFC,varid,'_FillValue',1.000000020040877e20);
    netcdf.putAtt(ncidFC,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
    netcdf.putAtt(ncidFC,varid,'institution','UCSD');
    netcdf.putAtt(ncidFC,varid,'wet weight:C ratio','9:1');

    netcdf.endDef(ncidFC);

    netcdf.putVar(ncidFC,vidlat,LAT);
    netcdf.putVar(ncidFC,vidlon,LON);
    netcdf.putVar(ncidFC,vidbioFC,AllFC);
    netcdf.putVar(ncidFC,vidtFC,yr);

    netcdf.close(ncidFC);

    %%
    ncdisp(file_tfc)

    %% tpc
    ncidPC = netcdf.create(file_tpc,cmode);

    time_dim = netcdf.defDim(ncidPC,'time',nt);
    lon_dim = netcdf.defDim(ncidPC,'nlon',ni);
    lat_dim = netcdf.defDim(ncidPC,'nlat',nj);

    vidtPC = netcdf.defVar(ncidPC,'time','NC_DOUBLE',time_dim);
    netcdf.putAtt(ncidPC,vidtPC,'long_name','time');
    netcdf.putAtt(ncidPC,vidtPC,'standard_name','time');
    netcdf.putAtt(ncidPC,vidtPC,'units','year' );
    netcdf.putAtt(ncidPC,vidtPC,'calendar','365_day');
    netcdf.putAtt(ncidPC,vidtPC,'axis','T');

    vidlon = netcdf.defVar(ncidPC,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
    netcdf.putAtt(ncidPC,vidlon,'long_name','longitude');
    netcdf.putAtt(ncidPC,vidlon,'standard_name','longitude');
    netcdf.putAtt(ncidPC,vidlon,'axis','X');

    vidlat = netcdf.defVar(ncidPC,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
    netcdf.putAtt(ncidPC,vidlat,'long_name','latitude');
    netcdf.putAtt(ncidPC,vidlat,'standard_name','latitude');
    netcdf.putAtt(ncidPC,vidlat,'axis','Y');

    vidbioPC = netcdf.defVar(ncidPC,'tpc','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
    netcdf.defVarChunking(ncidPC,vidbioPC,'CHUNKED',[10, 10, 1]);
    netcdf.putAtt(ncidPC,vidbioPC,'long_name','Total Annual Catch of Large Pelagic Fish');
    netcdf.putAtt(ncidPC,vidbioPC,'units','gWW' );
    netcdf.defVarFill(ncidPC,vidbioPC,false,1.000000020040877e20);
    netcdf.putAtt(ncidPC,vidbioPC,'missing value',1.000000020040877e20);

    varid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncidPC,varid,'creation_date',datestr(now));
    netcdf.putAtt(ncidPC,varid,'_FillValue',1.000000020040877e20);
    netcdf.putAtt(ncidPC,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
    netcdf.putAtt(ncidPC,varid,'institution','UCSD');
    netcdf.putAtt(ncidPC,varid,'wet weight:C ratio','9:1');

    netcdf.endDef(ncidPC);

    %tpb = single(tpb);
    netcdf.putVar(ncidPC,vidlat,LAT);
    netcdf.putVar(ncidPC,vidlon,LON);
    netcdf.putVar(ncidPC,vidbioPC,AllPC);
    netcdf.putVar(ncidPC,vidtPC,yr);

    netcdf.close(ncidPC);

    %% tdc
    ncidDC = netcdf.create(file_tdc,cmode);

    time_dim = netcdf.defDim(ncidDC,'time',nt);
    lon_dim = netcdf.defDim(ncidDC,'nlon',ni);
    lat_dim = netcdf.defDim(ncidDC,'nlat',nj);

    vidtDC = netcdf.defVar(ncidDC,'time','NC_DOUBLE',time_dim);
    netcdf.putAtt(ncidDC,vidtDC,'long_name','time');
    netcdf.putAtt(ncidDC,vidtDC,'standard_name','time');
    netcdf.putAtt(ncidDC,vidtDC,'calendar','365_day');
    netcdf.putAtt(ncidDC,vidtDC,'axis','T');
    netcdf.putAtt(ncidDC,vidtDC,'units','year' );

    vidlon = netcdf.defVar(ncidDC,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
    netcdf.putAtt(ncidDC,vidlon,'long_name','longitude');
    netcdf.putAtt(ncidDC,vidlon,'standard_name','longitude');
    netcdf.putAtt(ncidDC,vidlon,'axis','X');

    vidlat = netcdf.defVar(ncidDC,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
    netcdf.putAtt(ncidDC,vidlat,'long_name','latitude');
    netcdf.putAtt(ncidDC,vidlat,'standard_name','latitude');
    netcdf.putAtt(ncidDC,vidlat,'axis','Y');

    vidbioDC = netcdf.defVar(ncidDC,'tdc','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
    netcdf.defVarChunking(ncidDC,vidbioDC,'CHUNKED',[10, 10, 1]);
    netcdf.putAtt(ncidDC,vidbioDC,'long_name','Total Annual Catch of Demersal Fish');
    netcdf.putAtt(ncidDC,vidbioDC,'units','gWW' );
    netcdf.defVarFill(ncidDC,vidbioDC,false,1.000000020040877e20);
    netcdf.putAtt(ncidDC,vidbioDC,'missing value',1.000000020040877e20);

    varid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncidDC,varid,'creation_date',datestr(now));
    netcdf.putAtt(ncidDC,varid,'_FillValue',1.000000020040877e20);
    netcdf.putAtt(ncidDC,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
    netcdf.putAtt(ncidDC,varid,'institution','UCSD');
    netcdf.putAtt(ncidDC,varid,'wet weight:C ratio','9:1');

    netcdf.endDef(ncidDC);

    %tdb = single(tdb);
    netcdf.putVar(ncidDC,vidlat,LAT);
    netcdf.putVar(ncidDC,vidlon,LON);
    netcdf.putVar(ncidDC,vidbioDC,AllDC);
    netcdf.putVar(ncidDC,vidtDC,yr);

    netcdf.close(ncidDC);

end
