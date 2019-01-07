%%%%!! RUN SPINUP FOR ONE LOCATION
function TestonelocP()

global DAYS GRD NX
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a fcrit
global bent_eff rfrac CC
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel LPsel LDsel efn cfn

%fracm = 0.1:0.1:0.5;
Fmort = [0.0:0.1:1.0]; %[1.2:0.2:2.0]; %
RE = [1.0,0.5,0.1,0.05,0.01];
%BE = 0.05:0.05:0.1;
%CarCap = 0.5:0.5:2.0;
%encs = [1:4];
%cmaxs = [0;2;4];
%mets = 2:4;


% for E = 1:length(encs)
    efn = 2;%encs(E);
    
% for C = 1:length(cmaxs)
    cfn = 2;%cmaxs(C);
    
%     for M = 1:length(mets)
%         mfn = mets(M);
        
        for R = 1:length(RE)
            rfrac = RE(R);
            
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
                
                %! Setup spinup (loop last year of COBALT)
                load('/Volumes/GFDL/POEM_JLD/esm2m_hist/Data_ESM2Mhist_2000.mat');
                
                %! How long to run the model
                YEARS = 100;
                DAYS = 365;
                MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
                
                %! Where to run the model
                load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Data_grid_hindcast_NOTflipped.mat');
                ids = [40319,42639,41782,36334,38309,42744,30051,41284,38003,19327,20045];
                names = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
                
                %! Create a directory for output
                tfcrit = num2str(int64(100*fcrit));
                tld = num2str(1000+int64(100*LD_phi_MF));
                tkap = num2str(100+int64(10*K_a));
                tbe = num2str(100+int64(100*bent_eff));
                tmort = num2str(MORT);
                tcc = num2str(1000+int64(100*CC));
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
                    if (LPsel == 1 && LDsel == 1)
                        sel = 'L';
                    elseif (LPsel == 1)
                        sel = 'LP';
                    elseif (LDsel == 1)
                        sel = 'LD';
                    end
                end
                msel = num2str(100+int64(10*MPsel));
                if (pdc == 0)
                    coup = 'NoDc';
                elseif (pdc == 1)
                    coup = 'Dc';
                elseif (pdc == 2)
                    coup = 'PDc';
                end
                tcfn = num2str(cfn);
                tefn = num2str(efn);
                if (harv==1)
                    %simname = ['Ponly',coup,'_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit',tfcrit,'_D',tld(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_RE',tre(2:end),'_',sel,'_fish',tfish(2:end)];
                    %simname = ['Ponly',coup,'_TrefO_Hold_cmax-metab_MFeqMP_fcrit',tfcrit,'_D',tld(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_RE',tre(2:end),'_',sel,'_fish',tfish(2:end)];
                    simname = ['Ponly',coup,'_TrefO_cmax-metab',tcfn,'_enc',tefn,'_MFeqMP_fcrit',tfcrit,'_nmort',tmort,'_K',tkap(2:end),'_BE',tbe(2:end),'_CC',tcc(2:end),'_RE',tre(2:end),'_',sel,'_fish',tfish(2:end),'_Msel',msel(2:end)];
                else
                    %simname = ['Ponly',coup,'_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit',tfcrit,'_D',tld(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_RE',tre(2:end)];
                    %simname = ['Ponly',coup,'_TrefO_Hold_cmax-metab_MFeqMP_fcrit',tfcrit,'_D',tld(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_RE',tre(2:end)];
                    simname = ['Ponly',coup,'_TrefO_cmax-metab',tcfn,'_enc',tefn,'_MFeqMP_fcrit',tfcrit,'_nmort',tmort,'_K',tkap(2:end),'_BE',tbe(2:end),'_CC',tcc(2:end),'_RE',tre(2:end)];
                end
                if (~isdir(['/Volumes/GFDL/CSV/Matlab_new_size/',simname]))
                    mkdir(['/Volumes/GFDL/CSV/Matlab_new_size/',simname])
                end
%                 if (~isdir(['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_Big_sizes/',simname]))
%                     mkdir(['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_Big_sizes/',simname])
%                 end
                
                for L = 11;%1:length(ids);
                    ID = ids(L);
                    loc = names{L};
                    
                    NX = length(ID);
                    %! Initialize
                    [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,DAYS);
                    ENVR = sub_init_env(ID);
                    
                    %! Storage
                    Spinup_Sml_p = NaN*ones(DAYS,25);
                    Spinup_Med_p = NaN*ones(DAYS,25);
                    Spinup_Lrg_p = NaN*ones(DAYS,25);
                    Spinup_Cobalt = NaN*ones(DAYS,2);
                    
                    S_Sml_p = NaN*ones(12*YEARS,25);
                    S_Med_p = NaN*ones(12*YEARS,25);
                    S_Lrg_p = NaN*ones(12*YEARS,25);
                    S_Cobalt = NaN*ones(12*YEARS,2);
                    
                    %! Iterate forward in time with NO fishing
                    MNT=0;
                    for YR = 1:YEARS % years
                        for DAY = 1:DT:DAYS % days
                            
                            %%%! ticker
                            DY = int64(ceil(DAY));
                            [num2str(YR),' , ', num2str(mod(DY,365))]
                            
                            %%%! Future time step
                            [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
                                sub_futbioP(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,...
                                Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,dfrate,CC);
                            
                            %! Store last year
                            %if (YR==YEARS)
                                Spinup_Sml_p(DY,:) = [Sml_p.bio Sml_p.enc_f Sml_p.enc_p Sml_p.enc_d Sml_p.enc_zm Sml_p.enc_zl Sml_p.enc_be Sml_p.con_f Sml_p.con_p Sml_p.con_d Sml_p.con_zm Sml_p.con_zl Sml_p.con_be Sml_p.I Sml_p.nu Sml_p.gamma Sml_p.die Sml_p.rep Sml_p.rec Sml_p.clev Sml_p.prod Sml_p.pred Sml_p.nmort Sml_p.met Sml_p.caught];
                                Spinup_Med_p(DY,:) = [Med_p.bio Med_p.enc_f Med_p.enc_p Med_p.enc_d Med_p.enc_zm Med_p.enc_zl Med_p.enc_be Med_p.con_f Med_p.con_p Med_p.con_d Med_p.con_zm Med_p.con_zl Med_p.con_be Med_p.I Med_p.nu Med_p.gamma Med_p.die Med_p.rep Med_p.rec Med_p.clev Med_p.prod Med_p.pred Med_p.nmort Med_p.met Med_p.caught];
                                Spinup_Lrg_p(DY,:) = [Lrg_p.bio Lrg_p.enc_f Lrg_p.enc_p Lrg_p.enc_d Lrg_p.enc_zm Lrg_p.enc_zl Lrg_p.enc_be Lrg_p.con_f Lrg_p.con_p Lrg_p.con_d Lrg_p.con_zm Lrg_p.con_zl Lrg_p.con_be Lrg_p.I Lrg_p.nu Lrg_p.gamma Lrg_p.die Lrg_p.rep Lrg_p.rec Lrg_p.clev Lrg_p.prod Lrg_p.pred Lrg_p.nmort Lrg_p.met Lrg_p.caught];
                                Spinup_Cobalt(DY,:) = [ENVR.fZm ENVR.fZl];
                            %end
                            
                        end %Days
                        
                        %! Calculate monthly means and save
                        aa = (cumsum(MNTH)+1);
                        a = [1,aa(1:end-1)]; % start of the month
                        b = cumsum(MNTH); % end of the month
                        for i = 1:12
                            MNT = MNT+1; % Update monthly ticker
                            S_Cobalt(MNT,:) = mean(Spinup_Cobalt(a(i):b(i),:),1);
                            S_Sml_p(MNT,:) = mean(Spinup_Sml_p(a(i):b(i),:),1);
                            S_Med_p(MNT,:) = mean(Spinup_Med_p(a(i):b(i),:),1);
                            S_Lrg_p(MNT,:) = mean(Spinup_Lrg_p(a(i):b(i),:),1);
                        end
                        
                    end %Years
                    
                    %%% Save
                    csvwrite(['/Volumes/GFDL/CSV/Matlab_new_size/' simname '/Spinup_' loc '_Sml_p.csv'],S_Sml_p)
                    csvwrite(['/Volumes/GFDL/CSV/Matlab_new_size/' simname '/Spinup_' loc '_Med_p.csv'],S_Med_p)
                    csvwrite(['/Volumes/GFDL/CSV/Matlab_new_size/' simname '/Spinup_' loc '_Lrg_p.csv'],S_Lrg_p)
                    csvwrite(['/Volumes/GFDL/CSV/Matlab_new_size/' simname '/Spinup_' loc '_Cobalt.csv'],S_Cobalt)
                    
                end %Locations
            end %Fmort
        end %RE
    % end %met
% end %cmax
% end %enc
end
