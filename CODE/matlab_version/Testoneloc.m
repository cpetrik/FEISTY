%%%%!! RUN SPINUP FOR ONE LOCATION
function Testoneloc()

global Tref TrefP TrefB Dthresh SP DAYS GRD NX
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a fcrit
global bent_eff rfrac Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel LPsel LDsel

Fmort = 0.0; %[0.1:0.1:1.0];
CarCap = [0.25:0.25:3.0];
for C = 1%:length(CarCap)
    CC = CarCap(C);
    for F = 1:length(Fmort)
        %! Set fishing rate
        frate = Fmort(F);
        dfrate = Fmort(F)/365.0;
        %0=no fishing; 1=fishing
        if (frate>0)
            harv = 1;
        else
            harv = 0;
        end
        
        %! Make core parameters/constants (global)
        make_parameters()
        phen=0;
        
        %! Setup spinup (loop last year of COBALT)
        load('/Volumes/GFDL/POEM_JLD/esm2m_hist/Data_ESM2Mhist_2005.mat');
        
        %! Add phenology params from csv file with ID as row
        Tref = csvread('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/grid_phenol_T0raw_NOflip.csv'); %min temp for each yr at each location
        TrefP = Tref;
        TrefB = Tref;
        Dthresh = csvread('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/grid_phenol_DTraw_NOflip.csv');
        SP = csvread('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Gaussian_spawn_2mo.csv');
        
        %! How long to run the model
        YEARS = 100;
        DAYS = 365;
        
        %! Where to run the model
        load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Data_grid_hindcast_NOTflipped.mat');
        ids = [40319,42639,41782,36334,38309,42744,30051,41284,38003,19327,20045];
        names = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
        
        %! Create a directory for output
        tfcrit = num2str(int64(100*fcrit));
        tmz = num2str(100+int64(10*MF_phi_MZ));
        tld = num2str(1000+int64(100*LD_phi_MF));
        tbe = num2str(100+int64(100*bent_eff));
        tmort = num2str(MORT);
        tcc = num2str(100+int64(100*CC));
        if (rfrac >= 0.01)
            tre = num2str(10000+int64(1000*rfrac));
        else
            tre = num2str(100000+int64(round(10000*rfrac)));
        end
        if (frate >= 0.1)
            tfish = num2str(100+int64(10*frate));
        else
            tfish = num2str(1000+int64(100*frate));
        end
        if (MFsel == 1)
            if (LPsel == 1 && LDsel == 1)
                sel='All';
            else
                sel='MF';
            end
        else
            if (LPsel == 1)
                sel = 'LP';
            end
            if (LDsel == 1)
                sel = 'LD';
            end
        end
        if (pdc == 0)
            coup = 'NoDc';
        elseif (pdc == 1)
            coup = 'Dc';
        elseif (pdc == 2)
            coup = 'PDc';
        end
        if (harv==1)
            %simname = [coup,'_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit',tfcrit,'_BAassim','_nmort',tmort,'_BE',tbe(2:end),'_RE',tre(2:end),'_',sel,'_fish',tfish(2:end)];
            simname = [coup,'_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit',tfcrit,'_D',tld(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_RE',tre(2:end),'_',sel,'_fish',tfish(2:end)];
        else
            %simname = [coup,'_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit',tfcrit,'_BAassim','_nmort',tmort,'_BE',tbe(2:end),'_RE',tre(2:end)];
            simname = [coup,'_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit',tfcrit,'_D',tld(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_RE',tre(2:end)];
        end
        if (~isdir(['/Volumes/GFDL/CSV/Matlab_test_runs/',simname]))
            mkdir(['/Volumes/GFDL/CSV/Matlab_test_runs/',simname])
        end
        %     if (~isdir(['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/',simname]))
        %         mkdir(['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/',simname])
        %     end
        
        for L = 1%[1,2,6]
            ID = ids(L);
            loc = names{L};
            
            NX = length(ID);
            %! Initialize
            [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,phen,DAYS);
            Med_d.td(1) = 0.0;
            Lrg_d.td(1) = 0.0;
            ENVR = sub_init_env(ID);
            
            %! Storage
            Spinup_Sml_f = NaN*ones(YEARS*DAYS,28);
            Spinup_Sml_p = NaN*ones(YEARS*DAYS,28);
            Spinup_Sml_d = NaN*ones(YEARS*DAYS,28);
            Spinup_Med_f = NaN*ones(YEARS*DAYS,28);
            Spinup_Med_p = NaN*ones(YEARS*DAYS,28);
            Spinup_Med_d = NaN*ones(YEARS*DAYS,28);
            Spinup_Lrg_p = NaN*ones(YEARS*DAYS,28);
            Spinup_Lrg_d = NaN*ones(YEARS*DAYS,28);
            Spinup_Cobalt = NaN*ones(YEARS*DAYS,5);
            
            %! Iterate forward in time with NO fishing
            n=0;
            for YR = 1:YEARS % years
                %reset spawning flag
                if (phen == 1)
                    Med_f.S = zeros(NX,DAYS);
                    Lrg_d.S = zeros(NX,DAYS);
                    Lrg_p.S = zeros(NX,DAYS);
                end
                
                for DAY = 1:DT:DAYS % days
                    
                    %%%! ticker
                    n = n+1;
                    DY = int64(ceil(DAY));
                    [num2str(YR),' , ', num2str(mod(DY,365))]
                    
                    %%%! Future time step
                    [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
                        sub_futbio(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,...
                        Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,dfrate,CC);
                    
                    %! Store last year
                    %                 if (YR==YEARS)
                    %                     Spinup_Sml_f(DY,:) = [Sml_f.bio Sml_f.enc_f Sml_f.enc_p Sml_f.enc_d Sml_f.enc_zm ...
                    %                         Sml_f.enc_zl Sml_f.enc_be Sml_f.con_f Sml_f.con_p Sml_f.con_d Sml_f.con_zm ...
                    %                         Sml_f.con_zl Sml_f.con_be Sml_f.I Sml_f.nu Sml_f.gamma Sml_f.die Sml_f.rep ...
                    %                         Sml_f.rec Sml_f.egg Sml_f.clev Sml_f.DD Sml_f.S(DY) Sml_f.prod Sml_f.pred ...
                    %                         Sml_f.nmort Sml_f.met Sml_f.caught];
                    %                     Spinup_Sml_p(DY,:) = [Sml_p.bio Sml_p.enc_f Sml_p.enc_p Sml_p.enc_d Sml_p.enc_zm Sml_p.enc_zl Sml_p.enc_be Sml_p.con_f Sml_p.con_p Sml_p.con_d Sml_p.con_zm Sml_p.con_zl Sml_p.con_be Sml_p.I Sml_p.nu Sml_p.gamma Sml_p.die Sml_p.rep Sml_p.rec Sml_p.egg Sml_p.clev Sml_p.DD Sml_p.S(DY) Sml_p.prod Sml_p.pred Sml_p.nmort Sml_p.met Sml_p.caught];
                    %                     Spinup_Sml_d(DY,:) = [Sml_d.bio Sml_d.enc_f Sml_d.enc_p Sml_d.enc_d Sml_d.enc_zm Sml_d.enc_zl Sml_d.enc_be Sml_d.con_f Sml_d.con_p Sml_d.con_d Sml_d.con_zm Sml_d.con_zl Sml_d.con_be Sml_d.I Sml_d.nu Sml_d.gamma Sml_d.die Sml_d.rep Sml_d.rec Sml_d.egg Sml_d.clev Sml_d.DD Sml_d.S(DY) Sml_d.prod Sml_d.pred Sml_d.nmort Sml_d.met Sml_d.caught];
                    %                     Spinup_Med_f(DY,:) = [Med_f.bio Med_f.enc_f Med_f.enc_p Med_f.enc_d Med_f.enc_zm Med_f.enc_zl Med_f.enc_be Med_f.con_f Med_f.con_p Med_f.con_d Med_f.con_zm Med_f.con_zl Med_f.con_be Med_f.I Med_f.nu Med_f.gamma Med_f.die Med_f.rep Med_f.rec Med_f.egg Med_f.clev Med_f.DD Med_f.S(DY) Med_f.prod Med_f.pred Med_f.nmort Med_f.met Med_f.caught];
                    %                     Spinup_Med_p(DY,:) = [Med_p.bio Med_p.enc_f Med_p.enc_p Med_p.enc_d Med_p.enc_zm Med_p.enc_zl Med_p.enc_be Med_p.con_f Med_p.con_p Med_p.con_d Med_p.con_zm Med_p.con_zl Med_p.con_be Med_p.I Med_p.nu Med_p.gamma Med_p.die Med_p.rep Med_p.rec Med_p.egg Med_p.clev Med_p.DD Med_p.S(DY) Med_p.prod Med_p.pred Med_p.nmort Med_p.met Med_p.caught];
                    %                     Spinup_Med_d(DY,:) = [Med_d.bio Med_d.enc_f Med_d.enc_p Med_d.enc_d Med_d.enc_zm Med_d.enc_zl Med_d.enc_be Med_d.con_f Med_d.con_p Med_d.con_d Med_d.con_zm Med_d.con_zl Med_d.con_be Med_d.I Med_d.nu Med_d.gamma Med_d.die Med_d.rep Med_d.rec Med_d.egg Med_d.clev Med_d.DD Med_d.S(DY) Med_d.prod Med_d.pred Med_d.nmort Med_d.met Med_d.caught];
                    %                     Spinup_Lrg_p(DY,:) = [Lrg_p.bio Lrg_p.enc_f Lrg_p.enc_p Lrg_p.enc_d Lrg_p.enc_zm Lrg_p.enc_zl Lrg_p.enc_be Lrg_p.con_f Lrg_p.con_p Lrg_p.con_d Lrg_p.con_zm Lrg_p.con_zl Lrg_p.con_be Lrg_p.I Lrg_p.nu Lrg_p.gamma Lrg_p.die Lrg_p.rep Lrg_p.rec Lrg_p.egg Lrg_p.clev Lrg_p.DD Lrg_p.S(DY) Lrg_p.prod Lrg_p.pred Lrg_p.nmort Lrg_p.met Lrg_p.caught];
                    %                     Spinup_Lrg_d(DY,:) = [Lrg_d.bio Lrg_d.enc_f Lrg_d.enc_p Lrg_d.enc_d Lrg_d.enc_zm Lrg_d.enc_zl Lrg_d.enc_be Lrg_d.con_f Lrg_d.con_p Lrg_d.con_d Lrg_d.con_zm Lrg_d.con_zl Lrg_d.con_be Lrg_d.I Lrg_d.nu Lrg_d.gamma Lrg_d.die Lrg_d.rep Lrg_d.rec Lrg_d.egg Lrg_d.clev Lrg_d.DD Lrg_d.S(DY) Lrg_d.prod Lrg_d.pred Lrg_d.nmort Lrg_d.met Lrg_d.caught];
                    %                     Spinup_Cobalt(DY,:) = [BENT.mass ENVR.fZm ENVR.fZl ENVR.fB];
                    %                 end
                    Spinup_Sml_f(n,:) = [Sml_f.bio Sml_f.enc_f Sml_f.enc_p Sml_f.enc_d Sml_f.enc_zm ...
                        Sml_f.enc_zl Sml_f.enc_be Sml_f.con_f Sml_f.con_p Sml_f.con_d Sml_f.con_zm ...
                        Sml_f.con_zl Sml_f.con_be Sml_f.I Sml_f.nu Sml_f.gamma Sml_f.die Sml_f.rep ...
                        Sml_f.rec Sml_f.egg Sml_f.clev Sml_f.DD Sml_f.S(DY) Sml_f.prod Sml_f.pred ...
                        Sml_f.nmort Sml_f.met Sml_f.caught];
                    Spinup_Sml_p(n,:) = [Sml_p.bio Sml_p.enc_f Sml_p.enc_p Sml_p.enc_d Sml_p.enc_zm Sml_p.enc_zl Sml_p.enc_be Sml_p.con_f Sml_p.con_p Sml_p.con_d Sml_p.con_zm Sml_p.con_zl Sml_p.con_be Sml_p.I Sml_p.nu Sml_p.gamma Sml_p.die Sml_p.rep Sml_p.rec Sml_p.egg Sml_p.clev Sml_p.DD Sml_p.S(DY) Sml_p.prod Sml_p.pred Sml_p.nmort Sml_p.met Sml_p.caught];
                    Spinup_Sml_d(n,:) = [Sml_d.bio Sml_d.enc_f Sml_d.enc_p Sml_d.enc_d Sml_d.enc_zm Sml_d.enc_zl Sml_d.enc_be Sml_d.con_f Sml_d.con_p Sml_d.con_d Sml_d.con_zm Sml_d.con_zl Sml_d.con_be Sml_d.I Sml_d.nu Sml_d.gamma Sml_d.die Sml_d.rep Sml_d.rec Sml_d.egg Sml_d.clev Sml_d.DD Sml_d.S(DY) Sml_d.prod Sml_d.pred Sml_d.nmort Sml_d.met Sml_d.caught];
                    Spinup_Med_f(n,:) = [Med_f.bio Med_f.enc_f Med_f.enc_p Med_f.enc_d Med_f.enc_zm Med_f.enc_zl Med_f.enc_be Med_f.con_f Med_f.con_p Med_f.con_d Med_f.con_zm Med_f.con_zl Med_f.con_be Med_f.I Med_f.nu Med_f.gamma Med_f.die Med_f.rep Med_f.rec Med_f.egg Med_f.clev Med_f.DD Med_f.S(DY) Med_f.prod Med_f.pred Med_f.nmort Med_f.met Med_f.caught];
                    Spinup_Med_p(n,:) = [Med_p.bio Med_p.enc_f Med_p.enc_p Med_p.enc_d Med_p.enc_zm Med_p.enc_zl Med_p.enc_be Med_p.con_f Med_p.con_p Med_p.con_d Med_p.con_zm Med_p.con_zl Med_p.con_be Med_p.I Med_p.nu Med_p.gamma Med_p.die Med_p.rep Med_p.rec Med_p.egg Med_p.clev Med_p.DD Med_p.S(DY) Med_p.prod Med_p.pred Med_p.nmort Med_p.met Med_p.caught];
                    Spinup_Med_d(n,:) = [Med_d.bio Med_d.enc_f Med_d.enc_p Med_d.enc_d Med_d.enc_zm Med_d.enc_zl Med_d.enc_be Med_d.con_f Med_d.con_p Med_d.con_d Med_d.con_zm Med_d.con_zl Med_d.con_be Med_d.I Med_d.nu Med_d.gamma Med_d.die Med_d.rep Med_d.rec Med_d.egg Med_d.clev Med_d.DD Med_d.S(DY) Med_d.prod Med_d.pred Med_d.nmort Med_d.met Med_d.caught];
                    Spinup_Lrg_p(n,:) = [Lrg_p.bio Lrg_p.enc_f Lrg_p.enc_p Lrg_p.enc_d Lrg_p.enc_zm Lrg_p.enc_zl Lrg_p.enc_be Lrg_p.con_f Lrg_p.con_p Lrg_p.con_d Lrg_p.con_zm Lrg_p.con_zl Lrg_p.con_be Lrg_p.I Lrg_p.nu Lrg_p.gamma Lrg_p.die Lrg_p.rep Lrg_p.rec Lrg_p.egg Lrg_p.clev Lrg_p.DD Lrg_p.S(DY) Lrg_p.prod Lrg_p.pred Lrg_p.nmort Lrg_p.met Lrg_p.caught];
                    Spinup_Lrg_d(n,:) = [Lrg_d.bio Lrg_d.enc_f Lrg_d.enc_p Lrg_d.enc_d Lrg_d.enc_zm Lrg_d.enc_zl Lrg_d.enc_be Lrg_d.con_f Lrg_d.con_p Lrg_d.con_d Lrg_d.con_zm Lrg_d.con_zl Lrg_d.con_be Lrg_d.I Lrg_d.nu Lrg_d.gamma Lrg_d.die Lrg_d.rep Lrg_d.rec Lrg_d.egg Lrg_d.clev Lrg_d.DD Lrg_d.S(DY) Lrg_d.prod Lrg_d.pred Lrg_d.nmort Lrg_d.met Lrg_d.caught];
                    Spinup_Cobalt(n,:) = [BENT.mass BENT.pred ENVR.fZm ENVR.fZl ENVR.fB];
                end %Days
            end %Years
            
            %%% Save
            csvwrite(['/Volumes/GFDL/CSV/Matlab_test_runs/' simname '/Spinup_' loc '_Sml_f.csv'],Spinup_Sml_f)
            csvwrite(['/Volumes/GFDL/CSV/Matlab_test_runs/' simname '/Spinup_' loc '_Sml_p.csv'],Spinup_Sml_p)
            csvwrite(['/Volumes/GFDL/CSV/Matlab_test_runs/' simname '/Spinup_' loc '_Sml_d.csv'],Spinup_Sml_d)
            csvwrite(['/Volumes/GFDL/CSV/Matlab_test_runs/' simname '/Spinup_' loc '_Med_f.csv'],Spinup_Med_f)
            csvwrite(['/Volumes/GFDL/CSV/Matlab_test_runs/' simname '/Spinup_' loc '_Med_p.csv'],Spinup_Med_p)
            csvwrite(['/Volumes/GFDL/CSV/Matlab_test_runs/' simname '/Spinup_' loc '_Med_d.csv'],Spinup_Med_d)
            csvwrite(['/Volumes/GFDL/CSV/Matlab_test_runs/' simname '/Spinup_' loc '_Lrg_p.csv'],Spinup_Lrg_p)
            csvwrite(['/Volumes/GFDL/CSV/Matlab_test_runs/' simname '/Spinup_' loc '_Lrg_d.csv'],Spinup_Lrg_d)
            csvwrite(['/Volumes/GFDL/CSV/Matlab_test_runs/' simname '/Spinup_' loc '_Cobalt.csv'],Spinup_Cobalt)
            
        end %Locations
    end %Fmort
end %CC
end
