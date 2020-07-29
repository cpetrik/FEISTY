%%% Get CESM data
function ENVR = get_CESM(CESM,ID,DY)

    global GRD NX 

    %% Get data
    ENVR.Tp(:,1)  = CESM.Tp(ID,DY);
    ENVR.Tb(:,1)  = CESM.Tb(ID,DY);
    ENVR.Zm(:,1)  = CESM.Zm(ID,DY);
    ENVR.Zl(:,1)  = CESM.Zl(ID,DY);
    ENVR.det(:,1) = CESM.det(ID,DY);
%     ENVR.dZm(:,1) = CESM.dZm(ID,DY);
%     ENVR.dZl(:,1) = CESM.dZl(ID,DY);
%     ENVR.U(:,1)   = CESM.U(ID,DY);
%     ENVR.V(:,1)   = CESM.V(ID,DY);
    ENVR.fZm(:,1) = zeros(NX,1);
    ENVR.fZl(:,1) = zeros(NX,1);
    ENVR.fB(:,1)  = zeros(NX,1);
    ENVR.H(:,1)   = GRD.Z(ID);
%     ENVR.A(:,1)   = GRD.AREA(ID);
end
