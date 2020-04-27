%%%% File naming system
function [fname,simname] = sub_fname_ensemble6_samek()

global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac D J Sm A benc bcmx amet 
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MFsel MPsel MDsel LPsel LDsel Jsel efn cfn mfn
global tstep K CGRD ni nj frate dfrate kc ke

td = num2str(1000+int64(100*D));
tj = num2str(1000+int64(100*J));
tsm = num2str(1000+int64(100*Sm));
ta = num2str(1000+int64(100*A));
tbe = num2str(1000+int64(1000*bent_eff));
tmort = num2str(MORT);
tre = num2str(100000+int64(round(10000*rfrac)));
if (frate >= 0.1)
    tfish = num2str(1000+int64(100*frate));
    tF = num2str(1000+int64(100*frate*MFsel));
    tP = num2str(1000+int64(100*frate*LPsel));
    tD = num2str(1000+int64(100*frate*LDsel));
    tJ = num2str(1000+int64(100*Jsel));
else
    tfish = num2str(1000+int64(100*frate));
    tF = num2str(1000+int64(100*frate*MFsel));
    tP = num2str(1000+int64(100*frate*LPsel));
    tD = num2str(1000+int64(100*frate*LDsel));
    tJ = num2str(1000+int64(100*Jsel));
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
tmfn = num2str(1000+int64(100*amet));
tcfn = num2str(round(h));
tefn = num2str(round(gam));
tlam = num2str(1000+int64(1000*Lambda));
tkfn = num2str(1000+int64(1000*kt));
tke = num2str(1000+int64(1000*ke));
tkc = num2str(1000+int64(1000*kc));
tbfn = num2str(1000+int64(1000*bpow));
tbenc = num2str(1000+int64(1000*benc));
tbcmx = num2str(1000+int64(1000*bcmx));
tka = num2str(1000+int64(100*K_a));

simname = [coup,'_cmax',tcfn,'-b',tbcmx(2:end),'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_noCC_RE',tre(2:end),'_Ka',tka(2:end)];
ptext = ['Lam',tlam(2:end),'_enc',tefn,'-b',tbenc(2:end),'-k',tke(2:end),...
        '_cmax-k',tkc(2:end),'_met',tmfn(2:end),'-b',tbfn(2:end),'-k',tkfn(2:end)];

if (~isdir(['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',simname]))
    mkdir(['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',simname])
end


%! Setup netcdf path to store to
if (frate==0)
    fname = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',simname, '/Climatol_pristine_' ptext];
else
    fname  = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_' ptext];
end




end