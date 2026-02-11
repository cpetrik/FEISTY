%%% Get ESM data for global grid 
function ENVR = get_gESM(ESM,param,DY)

    %% Get data
    ENVR.Tp(:,1)  = ESM.Tp(param.ZID,DY);
    ENVR.Zm(:,1)  = ESM.Zm(param.ZID,DY);
    ENVR.Zl(:,1)  = ESM.Zl(param.ZID,DY);
    ENVR.dZm(:,1) = ESM.dZm(param.ZID,DY);
    ENVR.dZl(:,1) = ESM.dZl(param.ZID,DY);
    ENVR.Tb(:,1)  = ESM.Tb(DY);
    ENVR.det(:,1) = ESM.det(DY);
    ENVR.fZm(:,1) = zeros(param.NZ,1);
    ENVR.fZl(:,1) = zeros(param.NZ,1);
    ENVR.fB(:,1)  = 0.0;
    %ENVR.H(:,1)   = GRD.depth;
    ENVR.dz(:,1)  = ESM.dz(param.ZID,DY);
    %ENVR.A(:,1)   = GRD.area(param.ID);
end