%%% Get NEMURO data
function ENVR = get_NEMURO(NEMURO,ID,DY)

    global GRD NX 

    %% Get data
    ENVR.Tp(:,1)  = NEMURO.Tp(ID,DY);
    ENVR.Tb(:,1)  = NEMURO.Tb(ID,DY);
    ENVR.Zm(:,1)  = NEMURO.Zm(ID,DY);
    ENVR.Zl(:,1)  = NEMURO.Zl(ID,DY);
    ENVR.det(:,1) = NEMURO.det(ID,DY);
    ENVR.dZm(:,1) = NEMURO.dZm(ID,DY);
    ENVR.dZl(:,1) = NEMURO.dZl(ID,DY);
%     ENVR.U(:,1)   = NEMURO.U(ID,DY);
%     ENVR.V(:,1)   = NEMURO.V(ID,DY);
    ENVR.fZm(:,1) = zeros(NX,1);
    ENVR.fZl(:,1) = zeros(NX,1);
    ENVR.fB(:,1)  = zeros(NX,1);
    ENVR.H(:,1)   = GRD.Z(ID);
    ENVR.A(:,1)   = GRD.AREA(ID);
end
