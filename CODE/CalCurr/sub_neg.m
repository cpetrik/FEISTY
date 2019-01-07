%%% Negative value checks for interpolated COBALT values
function bio = sub_neg(bio)
    ID = (bio < 0);
    bio(ID) = 0.0;
    %ID2 = find(isnan(bio))
    %bio(ID2) = eps()
end
