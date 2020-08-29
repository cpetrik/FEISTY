%%% Get COBALT data
function ENVR = get_COBALT_1meso_empHPloss(COBALT,GRD,ID,DY)

    %% Get data
    ENVR.Tp(:,1)  = COBALT.Tp(ID,DY);
    ENVR.Tb(:,1)  = COBALT.Tb(ID,DY);
    ENVR.Zm(:,1)  = COBALT.Zm(ID,DY) + COBALT.Zl(ID,DY);
%     ENVR.dZm(:,1)  = COBALT.dZm(ID,DY) + COBALT.dZl(ID,DY);
    ENVR.det(:,1) = COBALT.det(ID,DY);
%     ENVR.U(:,1)   = COBALT.U(ID,DY);
%     ENVR.V(:,1)   = COBALT.V(ID,DY);
    ENVR.fZl(:,1) = zeros(length(ID),1);
    ENVR.fZb(:,1) = zeros(length(ID),1);
    ENVR.H(:,1)   = GRD.Z(ID);
    ENVR.A(:,1)   = GRD.AREA(ID);
    
    % HP loss is an empirical fitted fn of biomass and temp
    ENVR.dZm(:,1) = 10 .^ (-2.925 + 1.964.*log10(ENVR.Zm(:,1)+eps) + 1.958e-2.*ENVR.Tp(:,1));
    
end
