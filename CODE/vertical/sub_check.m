%%% Forward Euler checks
function bio = sub_check(bio)
    ID = (bio < 0);
    %ID = (bio < 0 | isnan(bio));
    bio(ID) = eps();

    ID2 = isnan(bio);
    bio(ID2) = 0.0;

end
