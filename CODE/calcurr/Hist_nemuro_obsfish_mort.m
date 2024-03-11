%%%%!! RUN HISTORIC WITH FISHING FOR ALL LOCATIONS
function Hist_nemuro_obsfish_mort()

vers = 'IPSL';

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
load('/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/IPSLdown/nemuro_ipsl_fmort_ID_annual_1980_2010_tempSc_assessment.mat',...
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
load('/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/IPSLdown/Data_grid_nemuro_ipsl.mat','GRD');
param.NX = length(GRD.Z);
param.ID = 1:param.NX;
NX = length(GRD.Z);
ID = 1:NX;

%! Create a directory for output
[fname,simname,outdir] = fname_hist(param,vers);

%! Storage variables
S_Sml_f_prod = zeros(NX,DAYS);
 S_Sml_p_prod = zeros(NX,DAYS);
 S_Sml_d_prod = zeros(NX,DAYS);
 S_Med_f_prod = zeros(NX,DAYS);
  S_Med_p_prod = zeros(NX,DAYS);
 S_Med_d_prod = zeros(NX,DAYS);
 S_Lrg_p_prod = zeros(NX,DAYS);
 S_Lrg_d_prod = zeros(NX,DAYS);

 S_Sml_f_die = zeros(NX,DAYS);
 S_Sml_p_die = zeros(NX,DAYS);
 S_Sml_d_die = zeros(NX,DAYS);
  S_Med_f_die = zeros(NX,DAYS);
 S_Med_p_die = zeros(NX,DAYS);
 S_Med_d_die = zeros(NX,DAYS);
 S_Lrg_p_die = zeros(NX,DAYS);
  S_Lrg_d_die = zeros(NX,DAYS);


%! Initialize
%!From a previous run
load([outdir 'Last_mo_Spinup_IPSL_All_fishobs_',simname,'.mat']);
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_hist(ID,DAYS,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_sml_f = [fname,'_mort_sml_f.nc'];
file_sml_p = [fname,'_mort_sml_p.nc'];
file_sml_d = [fname,'_mort_sml_d.nc'];
file_med_f = [fname,'_mort_med_f.nc'];
file_med_p = [fname,'_mort_med_p.nc'];
file_med_d = [fname,'_mort_med_d.nc'];
file_lrg_p = [fname,'_mort_lrg_p.nc'];
file_lrg_d = [fname,'_mort_lrg_d.nc'];

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
vidpredSF   = netcdf.defVar(ncidSF,'pred','double',[xy_dim,time_dim]);
viddieSF    = netcdf.defVar(ncidSF,'die','double',[xy_dim,time_dim]);
% vidclevSF   = netcdf.defVar(ncidSF,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidpredSP   = netcdf.defVar(ncidSP,'pred','double',[xy_dim,time_dim]);
viddieSP    = netcdf.defVar(ncidSP,'die','double',[xy_dim,time_dim]);
% vidclevSP   = netcdf.defVar(ncidSP,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidpredSD   = netcdf.defVar(ncidSD,'pred','double',[xy_dim,time_dim]);
viddieSD    = netcdf.defVar(ncidSD,'die','double',[xy_dim,time_dim]);
% vidclevSD   = netcdf.defVar(ncidSD,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidpredMF   = netcdf.defVar(ncidMF,'pred','double',[xy_dim,time_dim]);
viddieMF    = netcdf.defVar(ncidMF,'die','double',[xy_dim,time_dim]);
% vidclevMF   = netcdf.defVar(ncidMF,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidpredMP   = netcdf.defVar(ncidMP,'pred','double',[xy_dim,time_dim]);
viddieMP    = netcdf.defVar(ncidMP,'die','double',[xy_dim,time_dim]);
% vidclevMP   = netcdf.defVar(ncidMP,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidpredMD   = netcdf.defVar(ncidMD,'pred','double',[xy_dim,time_dim]);
viddieMD    = netcdf.defVar(ncidMD,'die','double',[xy_dim,time_dim]);
% vidclevMD   = netcdf.defVar(ncidMD,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidpredLP   = netcdf.defVar(ncidLP,'pred','double',[xy_dim,time_dim]);
viddieLP    = netcdf.defVar(ncidLP,'die','double',[xy_dim,time_dim]);
% vidclevLP   = netcdf.defVar(ncidLP,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidpredLD   = netcdf.defVar(ncidLD,'pred','double',[xy_dim,time_dim]);
viddieLD    = netcdf.defVar(ncidLD,'die','double',[xy_dim,time_dim]);
% vidclevLD   = netcdf.defVar(ncidLD,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLD);

%% %%%%%%%%%%%%%%%%%%%% Run the Model
%! Run model with no fishing
MNT=0;
for YR = 1:YEARS % years
    
    MY = num2str(modyrs(YR))
    load(['/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/IPSLdown/Data_nemuro_ipsl_',MY,'.mat'],'ESM');

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
        
	S_Sml_f_pred(:,DY) = Sml_f.pred;
                 S_Sml_p_pred(:,DY) = Sml_p.pred;
                 S_Sml_d_pred(:,DY) = Sml_d.pred;
                 S_Med_f_pred(:,DY) = Med_f.pred;
                 S_Med_p_pred(:,DY) = Med_p.pred;
                 S_Med_d_pred(:,DY) = Med_d.pred;
                 S_Lrg_p_pred(:,DY) = Lrg_p.pred;
                 S_Lrg_d_pred(:,DY) = Lrg_d.pred;
        
                 S_Sml_f_die(:,DY) = Sml_f.die;
                 S_Sml_p_die(:,DY) = Sml_p.die;
                 S_Sml_d_die(:,DY) = Sml_d.die;
                 S_Med_f_die(:,DY) = Med_f.die;
                 S_Med_p_die(:,DY) = Med_p.die;
                 S_Med_d_die(:,DY) = Med_d.die;
                 S_Lrg_p_die(:,DY) = Lrg_p.die;
                 S_Lrg_d_die(:,DY) = Lrg_d.die;
        
    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker
        
        %! Put vars of netcdf file
        netcdf.putVar(ncidSF,vidpredSF,[0 MNT-1],[NX 1],mean(S_Sml_f_pred(:,a(i):b(i)),2));
         			netcdf.putVar(ncidSP,vidpredSP,[0 MNT-1],[NX 1],mean(S_Sml_p_pred(:,a(i):b(i)),2));
         			netcdf.putVar(ncidSD,vidpredSD,[0 MNT-1],[NX 1],mean(S_Sml_d_pred(:,a(i):b(i)),2));
         			netcdf.putVar(ncidMF,vidpredMF,[0 MNT-1],[NX 1],mean(S_Med_f_pred(:,a(i):b(i)),2));
         			netcdf.putVar(ncidMP,vidpredMP,[0 MNT-1],[NX 1],mean(S_Med_p_pred(:,a(i):b(i)),2));
         			netcdf.putVar(ncidMD,vidpredMD,[0 MNT-1],[NX 1],mean(S_Med_d_pred(:,a(i):b(i)),2));
         			netcdf.putVar(ncidLP,vidpredLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_pred(:,a(i):b(i)),2));
         			netcdf.putVar(ncidLD,vidpredLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_pred(:,a(i):b(i)),2));
        
                   netcdf.putVar(ncidSF,viddieSF,[0 MNT-1],[NX 1],mean(S_Sml_f_die(:,a(i):b(i)),2));
         			netcdf.putVar(ncidSP,viddieSP,[0 MNT-1],[NX 1],mean(S_Sml_p_die(:,a(i):b(i)),2));
         			netcdf.putVar(ncidSD,viddieSD,[0 MNT-1],[NX 1],mean(S_Sml_d_die(:,a(i):b(i)),2));
         			netcdf.putVar(ncidMF,viddieMF,[0 MNT-1],[NX 1],mean(S_Med_f_die(:,a(i):b(i)),2));
         			netcdf.putVar(ncidMP,viddieMP,[0 MNT-1],[NX 1],mean(S_Med_p_die(:,a(i):b(i)),2));
         			netcdf.putVar(ncidMD,viddieMD,[0 MNT-1],[NX 1],mean(S_Med_d_die(:,a(i):b(i)),2));
         			netcdf.putVar(ncidLP,viddieLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_die(:,a(i):b(i)),2));
         			netcdf.putVar(ncidLD,viddieLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_die(:,a(i):b(i)),2));
					
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
