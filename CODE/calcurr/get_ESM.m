%%% Get NEMURO data
function ENVR = get_ESM(ESM,ID,DY)

    global GRD NX 

    %% Get data
    ENVR.Tp(:,1)  = ESM.Tp(ID,DY);
    ENVR.Tb(:,1)  = ESM.Tb(ID,DY);
    ENVR.Zm(:,1)  = ESM.Zm(ID,DY);
    ENVR.Zl(:,1)  = ESM.Zl(ID,DY);
    ENVR.det(:,1) = ESM.det(ID,DY);
    ENVR.dZm(:,1) = ESM.dZm(ID,DY);
    ENVR.dZl(:,1) = ESM.dZl(ID,DY);
    ENVR.fZm(:,1) = zeros(NX,1);
    ENVR.fZl(:,1) = zeros(NX,1);
    ENVR.fB(:,1)  = zeros(NX,1);
    ENVR.H(:,1)   = GRD.Z(ID);
    
    ENVR.dZm(:,1) = ESM.dZm(ID,DY);
    ENVR.dZl(:,1) = ESM.dZl(ID,DY);
    
end
