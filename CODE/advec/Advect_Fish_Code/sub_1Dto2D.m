%%%% Fish 1-D vector of ocean only cells to 2-D matrix of full grid
function bio2D = sub_1Dto2D(GRD,bio1D,param)

bio2D = nan*ones(param.ni,param.nj);
bio2D(GRD.ID) = bio1D;

end