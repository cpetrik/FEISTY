%%%%!! RUN HISTORIC WITH FISHING FOR ALL LOCATIONS
function Hist_nemuro_obsfish_prod()

esms = {'IPSL','GFDL','HAD'};
esms2 = {'ipsl','gfdl','hadley'};
for mod = 2:3
    vers = esms{mod};
    esm = esms2{mod};

    if mod==1
        gcpath = '/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/IPSLdown/';
    elseif mod==2
        gcpath = '/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/GFDLdown/';
    elseif mod==3
        gcpath = '/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/HADdown/';
    end
    

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
load([gcpath 'nemuro_',esm,'_fmort_ID_annual_1980_2010_tempSc_assessment.mat'],...
    'fmD','fmF','fmP');

% Set fishing rate as 1st year for fname
param.frate = nan;
param.frateF = fmF(:,1);
param.frateP = fmP(:,1);
param.frateD = fmD(:,1);
param.dfrateF = param.frateF/365.0;
param.dfrateP = param.frateP/365.0;
param.dfrateD = param.frateD/365.0;

%! Make core parameters/constants (global)
param = make_parameters(param);

%! How long to run the model
modyrs = 1980:2010;
YEARS = length(modyrs);
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Grid (choose where and when to run the model)
load([gcpath 'Data_grid_nemuro_',esm,'.mat'],'GRD');
param.NX = length(GRD.Z);
param.ID = 1:param.NX;
NX = length(GRD.Z);
ID = 1:NX;

%! Create a directory for output
[fname,simname,outdir] = fname_hist(param,vers);

%! Storage variables
S_Sml_f_nu = zeros(NX,DAYS);
S_Sml_p_nu = zeros(NX,DAYS);
S_Sml_d_nu = zeros(NX,DAYS);
S_Med_f_nu = zeros(NX,DAYS);
S_Med_p_nu = zeros(NX,DAYS);
S_Med_d_nu = zeros(NX,DAYS);
S_Lrg_p_nu = zeros(NX,DAYS);
S_Lrg_d_nu = zeros(NX,DAYS);

S_Sml_f_prod = zeros(NX,DAYS);
S_Sml_p_prod = zeros(NX,DAYS);
S_Sml_d_prod = zeros(NX,DAYS);
S_Med_f_prod = zeros(NX,DAYS);
S_Med_p_prod = zeros(NX,DAYS);
S_Med_d_prod = zeros(NX,DAYS);
S_Lrg_p_prod = zeros(NX,DAYS);
S_Lrg_d_prod = zeros(NX,DAYS);


%! Initialize
%!From a previous run
load([outdir 'Last_mo_Spinup_',vers,'_All_fishobs_',simname,'.mat']);
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_hist(ID,DAYS,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_sml_f = [fname,'_prod_sml_f.nc'];
file_sml_p = [fname,'_prod_sml_p.nc'];
file_sml_d = [fname,'_prod_sml_d.nc'];
file_med_f = [fname,'_prod_med_f.nc'];
file_med_p = [fname,'_prod_med_p.nc'];
file_med_d = [fname,'_prod_med_d.nc'];
file_lrg_p = [fname,'_prod_lrg_p.nc'];
file_lrg_d = [fname,'_prod_lrg_d.nc'];

ncidSF = netcdf.create(file_sml_f,'NC_WRITE');
ncidSP = netcdf.create(file_sml_p,'NC_WRITE');
ncidSD = netcdf.create(file_sml_d,'NC_WRITE');
ncidMF = netcdf.create(file_med_f,'NC_WRITE');
ncidMP = netcdf.create(file_med_p,'NC_WRITE');
ncidMD = netcdf.create(file_med_d,'NC_WRITE');
ncidLP = netcdf.create(file_lrg_p,'NC_WRITE');
ncidLD = netcdf.create(file_lrg_d,'NC_WRITE');

%! Dims of netcdf file
nt = 12*YEARS;
netcdf.setDefaultFormat('NC_FORMAT_64BIT');

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~5 minutes ... ']
xy_dim      = netcdf.defDim(ncidSF,'nid',NX);
time_dim    = netcdf.defDim(ncidSF,'ntime',nt);
vidprodSF   = netcdf.defVar(ncidSF,'prod','double',[xy_dim,time_dim]);
vidnuSF     = netcdf.defVar(ncidSF,'nu','double',[xy_dim,time_dim]);
% vidgammaSF  = netcdf.defVar(ncidSF,'gamma','double',[xy_dim,time_dim]);
% viddieSF    = netcdf.defVar(ncidSF,'die','double',[xy_dim,time_dim]);
% vidclevSF   = netcdf.defVar(ncidSF,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidprodSP   = netcdf.defVar(ncidSP,'prod','double',[xy_dim,time_dim]);
vidnuSP     = netcdf.defVar(ncidSP,'nu','double',[xy_dim,time_dim]);
% vidgammaSP  = netcdf.defVar(ncidSP,'gamma','double',[xy_dim,time_dim]);
% viddieSP    = netcdf.defVar(ncidSP,'die','double',[xy_dim,time_dim]);
% vidclevSP   = netcdf.defVar(ncidSP,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidprodSD   = netcdf.defVar(ncidSD,'prod','double',[xy_dim,time_dim]);
vidnuSD     = netcdf.defVar(ncidSD,'nu','double',[xy_dim,time_dim]);
% vidgammaSD  = netcdf.defVar(ncidSD,'gamma','double',[xy_dim,time_dim]);
% viddieSD    = netcdf.defVar(ncidSD,'die','double',[xy_dim,time_dim]);
% vidclevSD   = netcdf.defVar(ncidSD,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidprodMF   = netcdf.defVar(ncidMF,'prod','double',[xy_dim,time_dim]);
vidnuMF     = netcdf.defVar(ncidMF,'nu','double',[xy_dim,time_dim]);
% vidgammaMF  = netcdf.defVar(ncidMF,'gamma','double',[xy_dim,time_dim]);
% vidrepMF    = netcdf.defVar(ncidMF,'rep','double',[xy_dim,time_dim]);
% viddieMF    = netcdf.defVar(ncidMF,'die','double',[xy_dim,time_dim]);
% vidclevMF   = netcdf.defVar(ncidMF,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidprodMP   = netcdf.defVar(ncidMP,'prod','double',[xy_dim,time_dim]);
vidnuMP     = netcdf.defVar(ncidMP,'nu','double',[xy_dim,time_dim]);
% vidgammaMP  = netcdf.defVar(ncidMP,'gamma','double',[xy_dim,time_dim]);
% viddieMP    = netcdf.defVar(ncidMP,'die','double',[xy_dim,time_dim]);
% vidclevMP   = netcdf.defVar(ncidMP,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidprodMD   = netcdf.defVar(ncidMD,'prod','double',[xy_dim,time_dim]);
vidnuMD     = netcdf.defVar(ncidMD,'nu','double',[xy_dim,time_dim]);
% vidgammaMD  = netcdf.defVar(ncidMD,'gamma','double',[xy_dim,time_dim]);
% viddieMD    = netcdf.defVar(ncidMD,'die','double',[xy_dim,time_dim]);
% vidclevMD   = netcdf.defVar(ncidMD,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidprodLP   = netcdf.defVar(ncidLP,'prod','double',[xy_dim,time_dim]);
vidnuLP     = netcdf.defVar(ncidLP,'nu','double',[xy_dim,time_dim]);
% vidgammaLP  = netcdf.defVar(ncidLP,'gamma','double',[xy_dim,time_dim]);
% vidrepLP    = netcdf.defVar(ncidLP,'rep','double',[xy_dim,time_dim]);
% viddieLP    = netcdf.defVar(ncidLP,'die','double',[xy_dim,time_dim]);
% vidclevLP   = netcdf.defVar(ncidLP,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidprodLD   = netcdf.defVar(ncidLD,'prod','double',[xy_dim,time_dim]);
vidnuLD     = netcdf.defVar(ncidLD,'nu','double',[xy_dim,time_dim]);
% vidgammaLD  = netcdf.defVar(ncidLD,'gamma','double',[xy_dim,time_dim]);
% vidrepLD    = netcdf.defVar(ncidLD,'rep','double',[xy_dim,time_dim]);
% viddieLD    = netcdf.defVar(ncidLD,'die','double',[xy_dim,time_dim]);
% vidclevLD   = netcdf.defVar(ncidLD,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLD);


%% %%%%%%%%%%%%%%%%%%%% Run the Model
%! Run model with no fishing
MNT=0;
for YR = 1:YEARS % years
    
    MY = num2str(modyrs(YR))
    load([gcpath 'Data_nemuro_',esm,'_',MY,'.mat'],'ESM');

    param.frateF = fmF(:,YR);
    param.frateP = fmP(:,YR);
    param.frateD = fmD(:,YR);
    param.dfrateF = param.frateF/365.0;
    param.dfrateP = param.frateP/365.0;
    param.dfrateD = param.frateD/365.0;
    
    for DY = 1:DAYS % days
        
        
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio(DY,ESM,GRD,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);
        
        S_Sml_f_nu(:,DY) = Sml_f.nu;
                 S_Sml_p_nu(:,DY) = Sml_p.nu;
                 S_Sml_d_nu(:,DY) = Sml_d.nu;
                 S_Med_f_nu(:,DY) = Med_f.nu;
                 S_Med_p_nu(:,DY) = Med_p.nu;
                 S_Med_d_nu(:,DY) = Med_d.nu;
                 S_Lrg_p_nu(:,DY) = Lrg_p.nu;
                 S_Lrg_d_nu(:,DY) = Lrg_d.nu;
        
                 S_Sml_f_prod(:,DY) = Sml_f.prod;
                 S_Sml_p_prod(:,DY) = Sml_p.prod;
                 S_Sml_d_prod(:,DY) = Sml_d.prod;
                 S_Med_f_prod(:,DY) = Med_f.prod;
                 S_Med_p_prod(:,DY) = Med_p.prod;
                 S_Med_d_prod(:,DY) = Med_d.prod;
                 S_Lrg_p_prod(:,DY) = Lrg_p.prod;
                 S_Lrg_d_prod(:,DY) = Lrg_d.prod;
        
        
    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker
        
        %! Put vars of netcdf file
        netcdf.putVar(ncidSF,vidprodSF,[0 MNT-1],[NX 1],mean(S_Sml_f_prod(:,a(i):b(i)),2));
        			netcdf.putVar(ncidSP,vidprodSP,[0 MNT-1],[NX 1],mean(S_Sml_p_prod(:,a(i):b(i)),2));
         			netcdf.putVar(ncidSD,vidprodSD,[0 MNT-1],[NX 1],mean(S_Sml_d_prod(:,a(i):b(i)),2));
         			netcdf.putVar(ncidMF,vidprodMF,[0 MNT-1],[NX 1],mean(S_Med_f_prod(:,a(i):b(i)),2));
         			netcdf.putVar(ncidMP,vidprodMP,[0 MNT-1],[NX 1],mean(S_Med_p_prod(:,a(i):b(i)),2));
         			netcdf.putVar(ncidMD,vidprodMD,[0 MNT-1],[NX 1],mean(S_Med_d_prod(:,a(i):b(i)),2));
         			netcdf.putVar(ncidLP,vidprodLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_prod(:,a(i):b(i)),2));
         			netcdf.putVar(ncidLD,vidprodLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_prod(:,a(i):b(i)),2));
					
					netcdf.putVar(ncidSF,vidnuSF,[0 MNT-1],[NX 1],mean(S_Sml_f_nu(:,a(i):b(i)),2));
         			netcdf.putVar(ncidSP,vidnuSP,[0 MNT-1],[NX 1],mean(S_Sml_p_nu(:,a(i):b(i)),2));
         			netcdf.putVar(ncidSD,vidnuSD,[0 MNT-1],[NX 1],mean(S_Sml_d_nu(:,a(i):b(i)),2));
         			netcdf.putVar(ncidMF,vidnuMF,[0 MNT-1],[NX 1],mean(S_Med_f_nu(:,a(i):b(i)),2));
        			netcdf.putVar(ncidMP,vidnuMP,[0 MNT-1],[NX 1],mean(S_Med_p_nu(:,a(i):b(i)),2));
         			netcdf.putVar(ncidMD,vidnuMD,[0 MNT-1],[NX 1],mean(S_Med_d_nu(:,a(i):b(i)),2));
         			netcdf.putVar(ncidLP,vidnuLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_nu(:,a(i):b(i)),2));
         			netcdf.putVar(ncidLD,vidnuLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_nu(:,a(i):b(i)),2));
        
        
    end %Monthly mean
    
end %Years

%! Close save
netcdf.close(ncidSF);
netcdf.close(ncidSP);
netcdf.close(ncidSD);
netcdf.close(ncidMF);
netcdf.close(ncidMP);
netcdf.close(ncidMD);
netcdf.close(ncidLP);
netcdf.close(ncidLD);

end
