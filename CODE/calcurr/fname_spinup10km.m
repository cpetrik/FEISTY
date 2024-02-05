%%%% File naming system
function [fname,simname] = sub_fname_spinup10km(frate)
global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac D J Sm A benc bcmx amet 
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel Jsel efn cfn mfn
global tstep K CGRD ni nj

td = num2str(1000+int64(100*LD_phi_MP));
tj = num2str(1000+int64(100*MP_phi_S));
tsm = num2str(1000+int64(100*MF_phi_MZ));
ta = num2str(1000+int64(100*LP_phi_MF));
tbe = num2str(100+int64(100*bent_eff));
tmort = num2str(MORT);
tre = num2str(100000+int64(round(10000*rfrac)));
if (frate >= 0.1)
    tfish = num2str(100+int64(10*frate));
    tF = num2str(1000+int64(100*frate*MFsel));
    tP = num2str(1000+int64(100*frate*LPsel));
    tD = num2str(1000+int64(100*frate*LDsel));
    tJ = num2str(100+int64(10*Jsel));
else
    tfish = num2str(1000+int64(100*frate));
    tF = num2str(1000+int64(100*frate*MFsel));
    tP = num2str(1000+int64(100*frate*LPsel));
    tD = num2str(1000+int64(100*frate*LDsel));
end
if (MFsel > 0)
    if (LPsel > 0 && LDsel > 0)
        sel='All';
    else
        sel='F';
    end
else
    if (LPsel > 0 && LDsel > 0)
        sel = 'L';
    elseif (LPsel > 0)
        sel = 'P';
    elseif (LDsel > 0)
        sel = 'D';
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

if (~isdir(['/Volumes/GFDL/NC/Matlab_new_size/',simname,'/CalCurr']))
    mkdir(['/Volumes/GFDL/NC/Matlab_new_size/',simname,'/CalCurr'])
end

%! Setup netcdf path to store to
if (frate==0)
    fname = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/CalCurr/Spinup10km_pristine'];
elseif (Jsel~=0.1)
    fname = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/CalCurr/Spinup10km_', sel,'_fish',tfish(2:end),'_Juve',tJ(2:end)];
elseif (MFsel~=LPsel)
    fname = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/CalCurr/Spinup10km_fish_F',tF(2:end),'_P',tP(2:end),'_D',tD(2:end)];
else
    fname  = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/CalCurr/Spinup10km_', sel,'_fish',tfish(2:end)];  
end



end