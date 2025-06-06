%%%% File naming system
function [fname,simname] = sub_fname_core_spin_fished_obs_server()
global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac D J Sm A benc bcmx amet
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel Jsel efn cfn mfn
global dfrateF dfrateP dfrateD

td = num2str(1000+int64(100*LD_phi_MP));
tj = num2str(1000+int64(100*MP_phi_S));
tsm = num2str(1000+int64(100*MF_phi_MZ));
ta = num2str(1000+int64(100*LP_phi_MF));
tbe = num2str(100+int64(100*bent_eff));
tmort = num2str(MORT);
tre = num2str(100000+int64(round(10000*rfrac)));
tJ = num2str(100+int64(10*Jsel));
if (nanmean(dfrateF) == 0)
    tF = '0';
else
    selF = num2str(1000+int64(100*MFsel));
    tF = ['obs' selF(2:end)];
end
if (nanmean(dfrateP) == 0)
    tP = '0';
else
    selP = num2str(1000+int64(100*LPsel));
    tP = ['obs' selP(2:end)];
end
if (nanmean(dfrateD) == 0)
    tD = '0';
else
    selD = num2str(1000+int64(100*LDsel));  
    tD = ['obs' selD(2:end)];
end
if (nanmean(dfrateF) > 0)
    if (nanmean(dfrateP) > 0 && nanmean(dfrateD) > 0)
        sel='All';
        tfish='obs';
    end
end
if (pdc == 0)
    coup = 'NoDc';
elseif (pdc == 1)
    coup = 'Dc';
elseif (pdc == 2)
    coup = 'PDc';
end
tmfn = num2str(amet);
tcfn = num2str(h);
tefn = num2str(round(gam));
tkfn = num2str(1000+int64(1000*kt));
tbfn = num2str(1000+int64(1000*bpow));
tbenc = num2str(1000+int64(1000*benc));
tbcmx = num2str(1000+int64(1000*bcmx));

if (isnan(cfn))
    simname = [coup,'_enc',tefn,'-b',tbenc(2:end),'_m',tmfn,'-b',tbfn(2:end),'-k',tkfn(2:end),'_c',tcfn,'-b',tbcmx(2:end),'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_noCC_RE',tre(2:end)];
else
    simname = [coup,'_efn',num2str(efn),'_mfn',num2str(mfn),'_cfn',num2str(cfn),'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_noCC_RE',tre(2:end)];
end

if (~isfolder(['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/',simname,'/CORE']))
    mkdir(['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/',simname,'/CORE'])
end

%! Setup netcdf path to store to
if (Jsel~=0.1)
    fname = ['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/',simname, '/CORE/Core_spin_', sel,'_fish_',tfish,'_Juve',tJ(2:end)];
elseif (MFsel~=LPsel)
    fname = ['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/',simname, '/CORE/Core_spin_fish_F',tF,'_P',tP,'_D',tD];
else
    fname  = ['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/',simname, '/CORE/Core_spin_', sel,'_fish_',tfish];
end



end