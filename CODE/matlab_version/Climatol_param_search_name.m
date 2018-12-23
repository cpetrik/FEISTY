%%%%!! RUN Climatol FOR ALL LOCATIONS
function Climatol_param_search_name()

global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac D J Sm A benc bcmx amet 
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel Jsel efn cfn mfn
global tstep K CGRD ni nj frate kc ke

%! Setup Climatol (loop 5-year climatology of ESM2.6-COBALT)
load('/Volumes/GFDL/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_daily.mat');

%! How long to run the model
YEARS = 150;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! choose where and when to run the model
Pdrpbx = '/Users/cpetrik/Dropbox/';
load('/Volumes/GFDL/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_grid.mat');
NX = length(GRD.Z);
ID = 1:NX;
%2D Grid for advect-diff
load([Pdrpbx 'Princeton/POEM_2.0/CODE/Data/Hindcast_cgrid_cp2D.mat']);
[ni,nj] = size(CGRD.mask);


%%%%%%%%%%%%%%% Initialize Model Variables
%! Choose parameters from other models of my own combo
%1=Kiorboe&Hirst, 2=Hartvig, 3=mizer, 4=JC15, NA=mine
cfn=nan;
efn=nan;
mfn=nan;

paramt = cell(25,1);
simt = cell(25,1);

%! Make core parameters/constants (global)
frate = 0.3;
dfrate = frate/365.0;
make_parameters()

% PARAMETER CALIBRATION
ktemp = 0.0705:0.01:0.1105;
bmp = 0.15:0.025:0.25;
bees = 0.025:0.025:0.15;

z=14;
for k=4:length(ktemp)
    for j=1:length(bmp)   
        z=z+1;
        %bent_eff = bees(j); 
        bpow = bmp(j); 
        kt = ktemp(k);     

        tkfn = num2str(1000+int64(1000*kt));
        tbfn = num2str(1000+int64(1000*bpow));
        tbe = num2str(100+int64(100*bent_eff));
        ptext = ['_m-b',tbfn(2:end),'-k',tkfn(2:end)];
        %ptext = ['_m-k',tkfn(2:end),'_BE',tbe(2:end)];
      
%     %! Set globals that are function of other parameters
%     dfrate = frate/365.0;
%     MF_phi_MZ = Sm;
%     MP_phi_MZ = Sm*J;
%     MP_phi_LZ = J;
%     MP_phi_S  = J;
%     LP_phi_MF = 1.0*A;
%     LD_phi_MF = D*A;
%     LD_phi_MP = D;
    
    %! Create a directory for output
    fname = sub_fname_param_search(ptext);
    oname = sub_fname_nodir(frate);
        
    paramt{z} = ptext;
    simt{z} = oname;
    
    %! Storage variables
    % Dims
    nt = 12;
    
    S_Bent_bio = zeros(NX,DAYS);
    
    S_Sml_f = zeros(NX,DAYS);
    S_Sml_p = zeros(NX,DAYS);
    S_Sml_d = zeros(NX,DAYS);
    S_Med_f = zeros(NX,DAYS);
    S_Med_p = zeros(NX,DAYS);
    S_Med_d = zeros(NX,DAYS);
    S_Lrg_p = zeros(NX,DAYS);
    S_Lrg_d = zeros(NX,DAYS);
    
    S_Med_f_fish = zeros(NX,DAYS);
    S_Med_p_fish = zeros(NX,DAYS);
    S_Med_d_fish = zeros(NX,DAYS);
    S_Lrg_p_fish = zeros(NX,DAYS);
    S_Lrg_d_fish = zeros(NX,DAYS);
    
    Spinup_Sml_f.bio = NaN*ones(NX,nt);
    Spinup_Sml_p.bio = NaN*ones(NX,nt);
    Spinup_Sml_d.bio = NaN*ones(NX,nt);
    Spinup_Med_f.bio = NaN*ones(NX,nt);
    Spinup_Med_p.bio = NaN*ones(NX,nt);
    Spinup_Med_d.bio = NaN*ones(NX,nt);
    Spinup_Lrg_p.bio = NaN*ones(NX,nt);
    Spinup_Lrg_d.bio = NaN*ones(NX,nt);
    Spinup_Bent.bio = NaN*ones(NX,nt);
    
    Spinup_Med_f.yield = NaN*ones(NX,nt);
    Spinup_Med_p.yield = NaN*ones(NX,nt);
    Spinup_Med_d.yield = NaN*ones(NX,nt);
    Spinup_Lrg_p.yield = NaN*ones(NX,nt);
    Spinup_Lrg_d.yield = NaN*ones(NX,nt);
    
    
    %! Initialize
    [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,DAYS);
    Med_d.td(1:NX) = 0.0;
    Lrg_d.td(1:NX) = 0.0;
    ENVR = sub_init_env(ID);
    
    
    %% %%%%%%%%%%%%%%%%%%%% Run the Model
    %! Run model with no fishing
    MNT=0;
    for YR = 1:YEARS % years
        
        for DAY = 1:DT:DAYS % days
            
            %%%! Future time step
            DY = int64(ceil(DAY));
            [num2str(YR),' , ', num2str(mod(DY,365))]
            [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
                sub_futbio(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,...
                Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,dfrate);
            
            if (YR==YEARS)
                S_Bent_bio(:,DY) = BENT.mass;
                
                S_Sml_f(:,DY) = Sml_f.bio;
                S_Sml_p(:,DY) = Sml_p.bio;
                S_Sml_d(:,DY) = Sml_d.bio;
                S_Med_f(:,DY) = Med_f.bio;
                S_Med_p(:,DY) = Med_p.bio;
                S_Med_d(:,DY) = Med_d.bio;
                S_Lrg_p(:,DY) = Lrg_p.bio;
                S_Lrg_d(:,DY) = Lrg_d.bio;
                
                S_Med_f_fish(:,DY) = Med_f.caught;
                S_Med_p_fish(:,DY) = Med_p.caught;
                S_Med_d_fish(:,DY) = Med_d.caught;
                S_Lrg_p_fish(:,DY) = Lrg_p.caught;
                S_Lrg_d_fish(:,DY) = Lrg_d.caught;
            end
            
        end %Days
        
    end %Years
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        %! Put vars of netcdf file
        Spinup_Bent.bio(:,i) = mean(S_Bent_bio(:,a(i):b(i)),2);
        Spinup_Sml_f.bio(:,i) = mean(S_Sml_f(:,a(i):b(i)),2);
        Spinup_Sml_p.bio(:,i) = mean(S_Sml_p(:,a(i):b(i)),2);
        Spinup_Sml_d.bio(:,i) = mean(S_Sml_d(:,a(i):b(i)),2);
        Spinup_Med_f.bio(:,i) = mean(S_Med_f(:,a(i):b(i)),2);
        Spinup_Med_p.bio(:,i) = mean(S_Med_p(:,a(i):b(i)),2);
        Spinup_Med_d.bio(:,i) = mean(S_Med_d(:,a(i):b(i)),2);
        Spinup_Lrg_p.bio(:,i) = mean(S_Lrg_p(:,a(i):b(i)),2);
        Spinup_Lrg_d.bio(:,i) = mean(S_Lrg_d(:,a(i):b(i)),2);
        
        Spinup_Med_f.yield(:,i) = mean(S_Med_f_fish(:,a(i):b(i)),2);
        Spinup_Med_p.yield(:,i) = mean(S_Med_p_fish(:,a(i):b(i)),2);
        Spinup_Med_d.yield(:,i) = mean(S_Med_d_fish(:,a(i):b(i)),2);
        Spinup_Lrg_p.yield(:,i) = mean(S_Lrg_p_fish(:,a(i):b(i)),2);
        Spinup_Lrg_d.yield(:,i) = mean(S_Lrg_d_fish(:,a(i):b(i)),2);
        
    end %Monthly mean
    %%% Save
    save([fname,'.mat'],...
        'Spinup_Sml_f','Spinup_Sml_p','Spinup_Sml_d','Spinup_Med_f',...
        'Spinup_Med_p','Spinup_Med_d','Spinup_Lrg_p','Spinup_Lrg_d',...
        'Spinup_Bent')
    
    end
end

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,'/param_sens/'];
save([nfile 'Climatol_All_fish03_means_param_search.mat'],'paramt','simt',...
    '-append')

end

