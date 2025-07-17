%%%% Fish 1-D vector of ocean only cells to 2-D matrix of full grid
function bio1D = sub_2Dto1D(GRD,bio2D)

bio1D = bio2D(GRD.ID);

end