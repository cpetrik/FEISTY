%%% Get NEMURO data
function ENVR = get_ESM_2meso_empHP(ESM,ID,DY)

    global GRD NX 

    %% Get data
    ENVR.Tp(:,1)  = ESM.Tp(ID,DY);
    ENVR.Tb(:,1)  = ESM.Tb(ID,DY);
    ENVR.Zm(:,1)  = ESM.Zm(ID,DY);
    ENVR.Zl(:,1)  = ESM.Zl(ID,DY);
    ENVR.det(:,1) = ESM.det(ID,DY);
    
    ENVR.dZm(:,1) = 10 .^ (-2.617 + 1.989.*log10(ESM.Zm(ID,DY)+eps) + 1.732e-2.*ESM.Tp(ID,DY));
    ENVR.dZl(:,1) = 10 .^ (-2.954 + 2.228.*log10(ESM.Zl(ID,DY)+eps) + 2.556e-2.*ESM.Tp(ID,DY));

    ENVR.fZm(:,1) = zeros(NX,1);
    ENVR.fZl(:,1) = zeros(NX,1);
    ENVR.fB(:,1)  = zeros(NX,1);
    ENVR.H(:,1)   = GRD.Z(ID);

end
