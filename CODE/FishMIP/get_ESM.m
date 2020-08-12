%%% Get CESM data
function ENVR = get_ESM(ESM,GRD,param,DY)

    %% Get data
    ENVR.Tp(:,1)  = ESM.Tp(param.ID,DY);
    ENVR.Tb(:,1)  = ESM.Tb(param.ID,DY);
    ENVR.Zm(:,1)  = ESM.Zm(param.ID,DY);
    ENVR.det(:,1) = ESM.det(param.ID,DY);
%     ENVR.dZm(:,1) = ESM.dZm(param.ID,DY);
    ENVR.fZm(:,1) = zeros(param.NX,1);
    ENVR.fB(:,1)  = zeros(param.NX,1);
    ENVR.H(:,1)   = GRD.Z(param.ID);
%     ENVR.A(:,1)   = GRD.AREA(param.ID);
end
