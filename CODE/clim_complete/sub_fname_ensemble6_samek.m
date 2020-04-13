%%%% File naming system
function fname = sub_fname_ensemble6_samek(param)

td = num2str(1000+int64(100* param.D));
tj = num2str(1000+int64(100* param.J));
tsm = num2str(1000+int64(100* param.Sm));
ta = num2str(1000+int64(100* param.A));
tbe = num2str(1000+int64(1000* param.bent_eff));
tmort = num2str( param.MORT);
tre = num2str(100000+int64(round(10000* param.rfrac)));
if ( param.frate >= 0.1)
    tfish = num2str(1000+int64(100* param.frate));
    tF = num2str(1000+int64(100* param.frate* param.MFsel));
    tP = num2str(1000+int64(100* param.frate* param.LPsel));
    tD = num2str(1000+int64(100* param.frate* param.LDsel));
    tJ = num2str(1000+int64(100* param.Jsel));
else
    tfish = num2str(1000+int64(100* param.frate));
    tF = num2str(1000+int64(100* param.frate* param.MFsel));
    tP = num2str(1000+int64(100* param.frate* param.LPsel));
    tD = num2str(1000+int64(100* param.frate* param.LDsel));
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
tmfn = num2str(1000+int64(100* param.amet));
tcfn = num2str(round( param.h));
tefn = num2str(round( param.gam));
tlam = num2str(1000+int64(1000* param.Lambda));
tkfn = num2str(1000+int64(1000* param.kt));
tke = num2str(1000+int64(1000* param.ke));
tkc = num2str(1000+int64(1000* param.kc));
tbfn = num2str(1000+int64(1000* param.bpow));
tbenc = num2str(1000+int64(1000* param.benc));
tbcmx = num2str(1000+int64(1000* param.bcmx));
tka = num2str(1000+int64(100* param.K_a));

%spath = '/home/cpetrik/FEISTY_clim/Output/';
spath = '/Volumes/FEISTY/FEISTY_clim/Output/';

simname = [coup,'_cmax',tcfn,'-b',tbcmx(2:end),'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_noCC_RE',tre(2:end),'_Ka',tka(2:end)];
ptext = ['Lam',tlam(2:end),'_enc',tefn,'-b',tbenc(2:end),'-k',tke(2:end),...
        '_cmax-k',tkc(2:end),'_met',tmfn(2:end),'-b',tbfn(2:end),'-k',tkfn(2:end)];

if (~isdir([spath,simname]))
    mkdir([spath,simname])
end


%! Setup netcdf path to store to
if (param.frate==0)
    fname = [spath,simname, '/Climatol_pristine_' ptext];
else
    fname = [spath,simname, '/Climatol_', sel,'_fish',tfish(2:end),'_' ptext];
end




end