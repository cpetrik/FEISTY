%%%% File naming system
function [fname,simname] = sub_fname_spin_ipsl(param,vers)
frate = param.frate;

td = num2str(1000+int64(100 * param.LD_phi_MP));
ta = num2str(1000+int64(100 * param.LP_phi_MF));
tbe = num2str(100+int64(100 * param.bent_eff));
tcc = num2str(param.CC);
tmort = num2str(param.MORT);
tre = num2str(100000+int64(round(10000 * param.rfrac)));
if (nanmean(param.dfrateF) == 0)
    tF = '0';
else
    selF = num2str(1000+int64(100*param.MFsel));
    tF = ['obs' selF(2:end)];
end
if (nanmean(param.dfrateP) == 0)
    tP = '0';
else
    selP = num2str(1000+int64(100*param.LPsel));
    tP = ['obs' selP(2:end)];
end
if (nanmean(param.dfrateD) == 0)
    tD = '0';
else
    selD = num2str(1000+int64(100*param.LDsel));  
    tD = ['obs' selD(2:end)];
end
if (nanmean(param.dfrateF) > 0)
    if (nanmean(param.dfrateP) > 0 && nanmean(param.dfrateD) > 0)
        sel='All';
        tfish='obs';
    end
end
if (param.pdc == 0)
    coup = 'NoDc';
elseif (param.pdc == 1)
    coup = 'Dc';
elseif (param.pdc == 2)
    coup = 'PDc';
end
tmfn = num2str(int64(100 * param.amet)); %num2str(param.amet);
tcfn = num2str(param.h);
tefn = num2str(round(param.gam));
tkfn = num2str(1000+int64(1000 * param.kt));
tbfn = num2str(1000+int64(1000 * param.bpow));
tbenc = num2str(1000+int64(1000 * param.benc));
tbcmx = num2str(1000+int64(1000 * param.bcmx));
tlam = num2str(1000+int64(1000 * param.Lambda));

if (param.CC==0)
    simname = [coup,'_Lam',tlam(2:end),'_enc',tefn,'-b',tbenc(2:end),'_m',tmfn,'-b',tbfn(2:end),'-k',tkfn(2:end),'_c',tcfn,'-b',tbcmx(2:end),'_D',td(2:end),'_A',ta(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_noCC_RE',tre(2:end)];
else
    simname = [coup,'_Lam',tlam(2:end),'_enc',tefn,'-b',tbenc(2:end),'_m',tmfn,'-b',tbfn(2:end),'-k',tkfn(2:end),'_c',tcfn,'-b',tbcmx(2:end),'_D',td(2:end),'_A',ta(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc,'_RE',tre(2:end)];
end

outdir = ['/Volumes/petrik-lab/Feisty/NC/NEMURO/',simname,'/',vers,'/'];

if (~isfolder(outdir))
    mkdir(outdir)
end

%! Setup netcdf path to store to
if (frate==0)
    fname = [outdir,'/Spinup_',vers];
elseif (param.Jsel~=0.1)
    fname = [outdir,'/Spinup_',vers,'_', sel,'_fish',tfish(2:end),'_Juve',tJ(2:end)];
elseif (param.MFsel~=param.LPsel)
    fname = [outdir,'/Spinup_',vers,'_fish_F',tF(2:end),'_P',tP(2:end),'_D',tD(2:end)];
else
    fname  = [outdir,'/Spinup_',vers,sel,'_fish',tfish(2:end)];  
end



end