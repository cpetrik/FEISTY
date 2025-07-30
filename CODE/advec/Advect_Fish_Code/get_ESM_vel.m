%%% Get ESM data
function ENVR = get_ESM_vel(ESM,GRD,param,DY)

    %% Get data
    ENVR.Tp(:,1)  = ESM.Tp(param.ID,DY);
    ENVR.Tb(:,1)  = ESM.Tb(param.ID,DY);
    ENVR.Zm(:,1)  = ESM.Zm(param.ID,DY);
    ENVR.Zl(:,1)  = ESM.Zl(param.ID,DY);
    ENVR.det(:,1) = ESM.det(param.ID,DY);
    ENVR.dZm(:,1) = ESM.dZm(param.ID,DY);
    ENVR.dZl(:,1) = ESM.dZl(param.ID,DY);
    ENVR.fZm(:,1) = zeros(param.NX,1);
    ENVR.fZl(:,1) = zeros(param.NX,1);
    ENVR.fB(:,1)  = zeros(param.NX,1);
    ENVR.H(:,1)   = GRD.Z(param.ID);
    %ENVR.A(:,1)   = GRD.area(param.ID);
    ENVR.U   = ESM.U(:,:,DY);
    ENVR.V   = ESM.V(:,:,DY);
end