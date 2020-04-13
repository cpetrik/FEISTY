%%%% File naming system
function [fname,simname] = sub_fname(frate,param)

td = num2str(1000+int64(100 * param.LD_phi_MP));
tj = num2str(1000+int64(100 * param.MP_phi_S));
tsm = num2str(1000+int64(100 * param.MF_phi_MZ));
ta = num2str(1000+int64(100 * param.LP_phi_MF));
tbe = num2str(100+int64(100 * param.bent_eff));
tmort = num2str(param.MORT);
tre = num2str(100000+int64(round(10000 * param.rfrac)));
if (frate >= 0.1)
    tfish = num2str(100+int64(10*frate));
    tF = num2str(1000+int64(100*frate * param.MFsel));
    tP = num2str(1000+int64(100*frate * param.LPsel));
    tD = num2str(1000+int64(100*frate * param.LDsel));
    tJ = num2str(1000+int64(100 * param.Jsel));
else
    tfish = num2str(1000+int64(100*frate));
    tF = num2str(1000+int64(100*frate * param.MFsel));
    tP = num2str(1000+int64(100*frate *param.LPsel));
    tD = num2str(1000+int64(100*frate *param.LDsel));
    tJ = num2str(1000+int64(100* param.Jsel));
end
if (param.MFsel > 0)
    if (param.LPsel > 0 && param.LDsel > 0)
        sel='All';
    else
        sel='F';
    end
else
    if (param.LPsel > 0 && param.LDsel > 0)
        sel = 'L';
    elseif (param.LPsel > 0)
        sel = 'P';
    elseif (param.LDsel > 0)
        sel = 'D';
    end
end
if (param.pdc == 0)
    coup = 'NoDc';
elseif (param.pdc == 1)
    coup = 'Dc';
elseif (param.pdc == 2)
    coup = 'PDc';
end
tmfn = num2str(param.amet);
tcfn = num2str(param.h);
tefn = num2str(round(param.gam));
tkfn = num2str(1000+int64(1000 * param.kt));
tbfn = num2str(1000+int64(1000 * param.bpow));
tbenc = num2str(1000+int64(1000 * param.benc));
tbcmx = num2str(1000+int64(1000 * param.bcmx));

simname = [coup,'_enc',tefn,'-b',tbenc(2:end),'_m',tmfn,'-b',tbfn(2:end),'-k',tkfn(2:end),'_c',tcfn,'-b',tbcmx(2:end),'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_noCC_RE',tre(2:end)];

%! Setup netcdf path to store to
spath = '/home/cpetrik/FEISTY_clim/Output/';

if (~isdir([spath,simname]))
    mkdir([spath,simname])
end

if (frate==0)
    fname = [spath,simname, '/Climatol_pristine'];
elseif (param.Jsel ~= 0.1)
    fname = [spath,simname, '/Climatol_', sel,'_fish',tfish(2:end),'_Juve',tJ(2:end)];
elseif (param.MFsel ~= param.LPsel)
    fname = [spath,simname, '/Climatol_fish_F',tF(2:end),'_P',tP(2:end),'_D',tD(2:end)];
else
    fname = [spath,simname, '/Climatol_', sel,'_fish',tfish(2:end)];
end



end