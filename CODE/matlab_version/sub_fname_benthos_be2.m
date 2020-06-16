%%%% File naming system
function [fname,simname] = sub_fname_benthos_be2(frate)
global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac benc bcmx amet
global Nat_mrt MORT
global MD_phi_BE LD_phi_BE
global MDsel LDsel Jsel efn cfn mfn
global kc ke 

tbe = num2str(100+int64(100*bent_eff));
tmort = num2str(MORT);
tre = num2str(100000+int64(round(10000*rfrac)));
if (frate >= 0.1)
    tfish = num2str(100+int64(10*frate));
    tD = num2str(1000+int64(100*frate*LDsel));
    tJ = num2str(1000+int64(100*Jsel));
else
    tfish = num2str(1000+int64(100*frate));
    tD = num2str(1000+int64(100*frate*LDsel));
    tJ = num2str(1000+int64(100*Jsel));
end
if (LDsel > 0)
    if (MDsel > 0)
        sel='All';
    else
        sel='Adult';
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

%! Simname
if (isnan(cfn))
    simname = [coup,'_enc',tefn,'-b',tbenc(2:end),'_m',tmfn,'-b',tbfn(2:end),'-k',tkfn(2:end),'_c',tcfn,'-b',tbcmx(2:end),'_nmort',tmort,'_2B_BE',tbe(2:end),'_noCC_RE',tre(2:end)];
else
    simname = [coup,'_efn',num2str(efn),'_mfn',num2str(mfn),'_cfn',num2str(cfn),'_nmort',tmort,'_2B_BE',tbe(2:end),'_noCC_RE',tre(2:end)];
end

%! Folder
if (~isfolder(['/Volumes/FEISTY/NC/Matlab_new_size/',simname]))
    mkdir(['/Volumes/FEISTY/NC/Matlab_new_size/',simname])
end

%! Folder + simname
if (frate==0)
    fname = ['/Volumes/FEISTY/NC/Matlab_new_size/',simname, '/Climatol_pristine'];
elseif (Jsel~=0.1)
    fname = ['/Volumes/FEISTY/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_Juve',tJ(2:end)];
else
    fname  = ['/Volumes/FEISTY/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end)];  
end



end