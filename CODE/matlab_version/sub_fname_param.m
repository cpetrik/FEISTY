%%%% File naming system
function fname = sub_fname_param(ptext)
global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac D J Sm A benc bcmx amet
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel Jsel efn cfn mfn
global tstep K CGRD ni nj frate kc ke

orig = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

tfish = num2str(100+int64(10*frate));
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

simname = [orig '/param_sens'];

if (~isdir(['/Volumes/GFDL/NC/Matlab_new_size/',simname]))
    mkdir(['/Volumes/GFDL/NC/Matlab_new_size/',simname])
end


%! Setup netcdf path to store to
if (frate==0)
    fname = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_' ptext];
elseif (frate~=0.3)
    fname  = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish03_' ptext];
else
    fname  = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_' ptext];
end



end