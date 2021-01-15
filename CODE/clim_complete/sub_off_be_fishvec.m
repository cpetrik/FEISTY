%%% Offline coupling
function ENVR = sub_off_be_fishvec(fish,ENVR)
    
    %total biomass consumption of B by each fn type (g/m2)
    con = squeeze(fish.con(:,:,3)) .* fish.bio;
    %total biomass consumption of B by all (g/m2)
    tcon = sum(con,2);
    %fraction of Bs consumed 
    ENVR.fB = tcon ./ fish.bio(:,3);
end
