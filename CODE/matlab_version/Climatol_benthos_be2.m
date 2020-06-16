%%%%!! RUN Climatol FOR ALL LOCATIONS
function Climatol_benthos_be2()

global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac benc bcmx amet
global Nat_mrt MORT
global MD_phi_BE LD_phi_MD LD_phi_BE
global MDsel LDsel Jsel efn cfn mfn
global kc ke 

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
frate = 0.0; %Fish(F);
dfrate = frate/365.0;

%! Choose parameters from other models of my own combo
%1=Kiorboe&Hirst, 2=Hartvig, 3=mizer, 4=JC15, NA=mine
cfn=nan;
efn=nan;
mfn=nan;

%! Make core parameters/constants (global)
make_params_benthos_be2()

%! Setup Climatol (loop 5-year climatology of ESM2.6-COBALT)
load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_daily.mat','COBALT');

%! How long to run the model
YEARS = 150;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Grid; choose where and when to run the model
load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_grid.mat','GRD');
NX = length(GRD.Z);
ID = 1:NX;

%! Create a directory for output
fname = sub_fname_benthos_be2(frate);


%! Initialize
[Sml_d,Med_d,Lrg_d,BENT] = sub_init_benthos_be2(ID,DAYS);
ENVR = sub_init_env(ID);


%! Storage variables
S_Bent_Sbio = zeros(NX,DAYS);
S_Bent_Mbio = zeros(NX,DAYS);

S_Sml_d = zeros(NX,DAYS);
S_Med_d = zeros(NX,DAYS);
S_Lrg_d = zeros(NX,DAYS);

S_Sml_d_prod = zeros(NX,DAYS);
S_Med_d_prod = zeros(NX,DAYS);
S_Lrg_d_prod = zeros(NX,DAYS);

% Dims
nt = 12;

Spinup_Sml_d.bio = NaN*ones(NX,nt);
Spinup_Med_d.bio = NaN*ones(NX,nt);
Spinup_Lrg_d.bio = NaN*ones(NX,nt);
Spinup_SBent.bio = NaN*ones(NX,nt);
Spinup_MBent.bio = NaN*ones(NX,nt);

%% %%%%%%%%%%%%%%%%%%%% Run the Model
%! Run model with no fishing
MNT=0;
for YR = 1:YEARS % years
    
    for DAY = 1:DT:DAYS % days
        
        %%%! Future time step
        DY = int64(ceil(DAY));
        [num2str(YR),' , ', num2str(mod(DY,365))]
        [Sml_d,Med_d,Lrg_d,BENT,ENVR] = sub_tbio_benthos_be2(ID,DY,COBALT,ENVR,...
            Sml_d,Med_d,Lrg_d,BENT,dfrate);
        
        if (YR==YEARS)
            S_Bent_Sbio(:,DY) = BENT.sm;
            S_Bent_Mbio(:,DY) = BENT.md;
            
            S_Sml_d(:,DY) = Sml_d.bio;
            S_Med_d(:,DY) = Med_d.bio;
            S_Lrg_d(:,DY) = Lrg_d.bio;
            
            S_Sml_d_prod(:,DY) = Sml_d.prod;
            S_Med_d_prod(:,DY) = Med_d.prod;
            S_Lrg_d_prod(:,DY) = Lrg_d.prod;
        end
        
    end %Days
    
end %Years

%! Calculate monthly means and save
aa = (cumsum(MNTH)+1);
a = [1,aa(1:end-1)]; % start of the month
b = cumsum(MNTH); % end of the month
for i = 1:12
    %! Put vars of netcdf file
    Spinup_SBent.bio(:,i) = mean(S_Bent_Sbio(:,a(i):b(i)),2);
    Spinup_MBent.bio(:,i) = mean(S_Bent_Mbio(:,a(i):b(i)),2);
    
    Spinup_Sml_d.bio(:,i) = mean(S_Sml_d(:,a(i):b(i)),2);
    Spinup_Med_d.bio(:,i) = mean(S_Med_d(:,a(i):b(i)),2);
    Spinup_Lrg_d.bio(:,i) = mean(S_Lrg_d(:,a(i):b(i)),2);
    
    Spinup_Sml_d.prod(:,i) = mean(S_Med_d_prod(:,a(i):b(i)),2);
    Spinup_Med_d.prod(:,i) = mean(S_Med_d_prod(:,a(i):b(i)),2);
    Spinup_Lrg_d.prod(:,i) = mean(S_Lrg_d_prod(:,a(i):b(i)),2);
    
end %Monthly mean

%%% Save
save([fname,'.mat'],...
    'Spinup_Sml_d','Spinup_Med_d','Spinup_Lrg_d',...
    'Spinup_SBent','Spinup_MBent');

end
