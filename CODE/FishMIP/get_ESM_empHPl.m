%%% Get ESM data
function ENVR = get_ESM_empHPl(ESM,GRD,param,DY)

    %% Get data
    ENVR.Tp(:,1)  = ESM.Tp(param.ID,DY);
    ENVR.Tb(:,1)  = ESM.Tb(param.ID,DY);
    ENVR.Zm(:,1)  = ESM.Zm(param.ID,DY);
    ENVR.det(:,1) = ESM.det(param.ID,DY);
    ENVR.fZl(:,1) = zeros(param.NX,1);
    ENVR.fB(:,1)  = zeros(param.NX,1);
    ENVR.H(:,1)   = GRD.Z(param.ID);

% HP loss is an empirical fitted fn of biomass and temp
    ENVR.dZm(:,1) = 10 .^ (-2.925 + 1.964.*log10(ENVR.Zm(:,1)+eps) + 1.958e-2.*ENVR.Tp(:,1));
    
end
