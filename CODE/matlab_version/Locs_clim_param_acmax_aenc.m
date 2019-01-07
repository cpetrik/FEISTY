%%%%!! RUN SPINUP FOR ONE LOCATION
function Locs_clim_param_acmax_aenc_search()

global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac D J Sm A benc bcmx amet
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel Jsel efn cfn mfn
global tstep K CGRD ni nj

%! Setup Climatol (loop 5-year climatology of ESM2.6-COBALT)
load('/Volumes/GFDL/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_daily.mat');

%! How long to run the model
YEARS = 150;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! choose where and when to run the model
load('/Volumes/GFDL/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_grid.mat');
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
names = abbrev;
ID = ids;

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
frate = 0.3; %Fish(F);
dfrate = frate/365.0;

%! Choose parameters from other models of my own combo
%1=Kiorboe&Hirst, 2=Hartvig, 3=mizer, 4=JC15, NA=mine
cfn=nan;
efn=nan;
mfn=nan;

%! Make core parameters/constants (global)
make_parameters_mid()
A = 0.5;
bpow = 0.175;
bcmx = 0.25;

aep = 10:10:100;
acp = 10:10:100;

for k=1:length(acp)
    for j=1:length(aep)   
        h = acp(k);         % coeff on Cmax
        gam = aep(j);       % coeff on search area
        
        %! Set globals that are function of other parameters
        dfrate = frate/365.0;
        MF_phi_MZ = Sm;
        MP_phi_MZ = Sm*J;
        MP_phi_LZ = J;
        MP_phi_S  = J;
        LP_phi_MF = 1.0*A;
        LD_phi_MF = D*A;
        LD_phi_MP = D;
        
        %! Create a directory for output
        fname = sub_fname_param_locs();
        
        NX = length(ID);
        %! Initialize
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,DAYS);
        Med_d.td(1) = 0.0;
        Lrg_d.td(1) = 0.0;
        ENVR = sub_init_env(ID);
        
        %! Storage
        Spinup_Sml_f = NaN*ones(DAYS,26,NX);
        Spinup_Sml_p = NaN*ones(DAYS,26,NX);
        Spinup_Sml_d = NaN*ones(DAYS,26,NX);
        Spinup_Med_f = NaN*ones(DAYS,26,NX);
        Spinup_Med_p = NaN*ones(DAYS,26,NX);
        Spinup_Med_d = NaN*ones(DAYS,26,NX);
        Spinup_Lrg_p = NaN*ones(DAYS,26,NX);
        Spinup_Lrg_d = NaN*ones(DAYS,26,NX);
        Spinup_Cobalt = NaN*ones(DAYS,5,NX);
        
        S_Sml_f = NaN*ones(12*YEARS,26,NX);
        S_Sml_p = NaN*ones(12*YEARS,26,NX);
        S_Sml_d = NaN*ones(12*YEARS,26,NX);
        S_Med_f = NaN*ones(12*YEARS,26,NX);
        S_Med_p = NaN*ones(12*YEARS,26,NX);
        S_Med_d = NaN*ones(12*YEARS,26,NX);
        S_Lrg_p = NaN*ones(12*YEARS,26,NX);
        S_Lrg_d = NaN*ones(12*YEARS,26,NX);
        S_Cobalt = NaN*ones(12*YEARS,5,NX);
        
        %! Iterate forward in time with NO fishing
        MNT=0;
        for YR = 1:YEARS % years
            num2str(YR)
            for DAY = 1:DT:DAYS % days
                
                %%%! ticker
                DY = int64(ceil(DAY));
                %[num2str(YR),' , ', num2str(mod(DY,365))]
                
                %%%! Future time step
                [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
                    sub_futbio(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,...
                    Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,dfrate);
                
                %! Store last year
                %if (YR==YEARS)
                Spinup_Sml_f(DY,:,:) = [Sml_f.bio Sml_f.enc_f Sml_f.enc_p Sml_f.enc_d Sml_f.enc_zm ...
                    Sml_f.enc_zl Sml_f.enc_be Sml_f.con_f Sml_f.con_p Sml_f.con_d Sml_f.con_zm ...
                    Sml_f.con_zl Sml_f.con_be Sml_f.I Sml_f.nu Sml_f.gamma Sml_f.die Sml_f.rep ...
                    Sml_f.rec Sml_f.clev Sml_f.prod Sml_f.pred Sml_f.nmort Sml_f.met Sml_f.caught Sml_f.fmort]';
                Spinup_Sml_p(DY,:,:) = [Sml_p.bio Sml_p.enc_f Sml_p.enc_p Sml_p.enc_d Sml_p.enc_zm Sml_p.enc_zl Sml_p.enc_be Sml_p.con_f Sml_p.con_p Sml_p.con_d Sml_p.con_zm Sml_p.con_zl Sml_p.con_be Sml_p.I Sml_p.nu Sml_p.gamma Sml_p.die Sml_p.rep Sml_p.rec Sml_p.clev Sml_p.prod Sml_p.pred Sml_p.nmort Sml_p.met Sml_p.caught Sml_p.fmort]';
                Spinup_Sml_d(DY,:,:) = [Sml_d.bio Sml_d.enc_f Sml_d.enc_p Sml_d.enc_d Sml_d.enc_zm Sml_d.enc_zl Sml_d.enc_be Sml_d.con_f Sml_d.con_p Sml_d.con_d Sml_d.con_zm Sml_d.con_zl Sml_d.con_be Sml_d.I Sml_d.nu Sml_d.gamma Sml_d.die Sml_d.rep Sml_d.rec Sml_d.clev Sml_d.prod Sml_d.pred Sml_d.nmort Sml_d.met Sml_d.caught Sml_d.fmort]';
                Spinup_Med_f(DY,:,:) = [Med_f.bio Med_f.enc_f Med_f.enc_p Med_f.enc_d Med_f.enc_zm Med_f.enc_zl Med_f.enc_be Med_f.con_f Med_f.con_p Med_f.con_d Med_f.con_zm Med_f.con_zl Med_f.con_be Med_f.I Med_f.nu Med_f.gamma Med_f.die Med_f.rep Med_f.rec Med_f.clev Med_f.prod Med_f.pred Med_f.nmort Med_f.met Med_f.caught Med_f.fmort]';
                Spinup_Med_p(DY,:,:) = [Med_p.bio Med_p.enc_f Med_p.enc_p Med_p.enc_d Med_p.enc_zm Med_p.enc_zl Med_p.enc_be Med_p.con_f Med_p.con_p Med_p.con_d Med_p.con_zm Med_p.con_zl Med_p.con_be Med_p.I Med_p.nu Med_p.gamma Med_p.die Med_p.rep Med_p.rec Med_p.clev Med_p.prod Med_p.pred Med_p.nmort Med_p.met Med_p.caught Med_p.fmort]';
                Spinup_Med_d(DY,:,:) = [Med_d.bio Med_d.enc_f Med_d.enc_p Med_d.enc_d Med_d.enc_zm Med_d.enc_zl Med_d.enc_be Med_d.con_f Med_d.con_p Med_d.con_d Med_d.con_zm Med_d.con_zl Med_d.con_be Med_d.I Med_d.nu Med_d.gamma Med_d.die Med_d.rep Med_d.rec Med_d.clev Med_d.prod Med_d.pred Med_d.nmort Med_d.met Med_d.caught Med_d.fmort]';
                Spinup_Lrg_p(DY,:,:) = [Lrg_p.bio Lrg_p.enc_f Lrg_p.enc_p Lrg_p.enc_d Lrg_p.enc_zm Lrg_p.enc_zl Lrg_p.enc_be Lrg_p.con_f Lrg_p.con_p Lrg_p.con_d Lrg_p.con_zm Lrg_p.con_zl Lrg_p.con_be Lrg_p.I Lrg_p.nu Lrg_p.gamma Lrg_p.die Lrg_p.rep Lrg_p.rec Lrg_p.clev Lrg_p.prod Lrg_p.pred Lrg_p.nmort Lrg_p.met Lrg_p.caught Lrg_p.fmort]';
                Spinup_Lrg_d(DY,:,:) = [Lrg_d.bio Lrg_d.enc_f Lrg_d.enc_p Lrg_d.enc_d Lrg_d.enc_zm Lrg_d.enc_zl Lrg_d.enc_be Lrg_d.con_f Lrg_d.con_p Lrg_d.con_d Lrg_d.con_zm Lrg_d.con_zl Lrg_d.con_be Lrg_d.I Lrg_d.nu Lrg_d.gamma Lrg_d.die Lrg_d.rep Lrg_d.rec Lrg_d.clev Lrg_d.prod Lrg_d.pred Lrg_d.nmort Lrg_d.met Lrg_d.caught Lrg_d.fmort]';
                Spinup_Cobalt(DY,:,:) = [BENT.mass BENT.pred ENVR.fZm ENVR.fZl ENVR.fB]';
                %Spinup_Cobalt(DY,:,:) = [BENT.sm BENT.md ENVR.fZm ENVR.fZl ENVR.fB]';
                %end
                
                
            end %Days
            
            %! Calculate monthly means and save
            aa = (cumsum(MNTH)+1);
            a = [1,aa(1:end-1)]; % start of the month
            b = cumsum(MNTH); % end of the month
            for i = 1:12
                MNT = MNT+1; % Update monthly ticker
                S_Cobalt(MNT,:,:) = mean(Spinup_Cobalt(a(i):b(i),:,:),1);
                S_Sml_f(MNT,:,:) = mean(Spinup_Sml_f(a(i):b(i),:,:),1);
                S_Sml_p(MNT,:,:) = mean(Spinup_Sml_p(a(i):b(i),:,:),1);
                S_Sml_d(MNT,:,:) = mean(Spinup_Sml_d(a(i):b(i),:,:),1);
                S_Med_f(MNT,:,:) = mean(Spinup_Med_f(a(i):b(i),:,:),1);
                S_Med_p(MNT,:,:) = mean(Spinup_Med_p(a(i):b(i),:,:),1);
                S_Med_d(MNT,:,:) = mean(Spinup_Med_d(a(i):b(i),:,:),1);
                S_Lrg_p(MNT,:,:) = mean(Spinup_Lrg_p(a(i):b(i),:,:),1);
                S_Lrg_d(MNT,:,:) = mean(Spinup_Lrg_d(a(i):b(i),:,:),1);
            end
            
        end %Years
        
        %%% Save
        save([fname '_locs.mat'],...
            'S_Sml_f','S_Sml_p','S_Sml_d','S_Med_f','S_Med_p','S_Med_d',...
            'S_Lrg_p','S_Lrg_d','S_Cobalt')
        
    end
end

end

