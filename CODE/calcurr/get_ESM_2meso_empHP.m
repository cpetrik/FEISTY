%%% Get NEMURO data
function ENVR = get_ESM_2meso_empHP(ESM,GRD,param,DY)

    
    %% Get data
    ENVR.Tp(:,1)  = ESM.Tp(param.ID,DY);
    ENVR.Tb(:,1)  = ESM.Tb(param.ID,DY);
    ENVR.Zm(:,1)  = ESM.Zm(param.ID,DY);
    ENVR.Zl(:,1)  = ESM.Zl(param.ID,DY);
    ENVR.det(:,1) = ESM.det(param.ID,DY);
    
    ENVR.dZm(:,1) = 10 .^ (-2.617 + 1.989.*log10(ESM.Zm(param.ID,DY)+eps) + 1.732e-2.*ESM.Tp(param.ID,DY));
    ENVR.dZl(:,1) = 10 .^ (-2.954 + 2.228.*log10(ESM.Zl(param.ID,DY)+eps) + 2.556e-2.*ESM.Tp(param.ID,DY));

    ENVR.fZm(:,1) = zeros(param.NX,1);
    ENVR.fZl(:,1) = zeros(param.NX,1);
    ENVR.fB(:,1)  = zeros(param.NX,1);
    ENVR.H(:,1)   = GRD.Z(param.ID);

end
