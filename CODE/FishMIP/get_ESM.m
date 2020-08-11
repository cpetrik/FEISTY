%%% Get CESM data
function ENVR = get_ESM(ESM,GRD,ID,DY)

    %% Get data
    ENVR.Tp(:,1)  = ESM.Tp(ID,DY);
    ENVR.Tb(:,1)  = ESM.Tb(ID,DY);
    ENVR.Zm(:,1)  = ESM.Zm(ID,DY);
    ENVR.det(:,1) = ESM.det(ID,DY);
%     ENVR.dZm(:,1) = ESM.dZm(ID,DY);
    ENVR.fZm(:,1) = zeros(NX,1);
    ENVR.fB(:,1)  = zeros(NX,1);
    ENVR.H(:,1)   = GRD.Z(ID);
%     ENVR.A(:,1)   = GRD.AREA(ID);
end
