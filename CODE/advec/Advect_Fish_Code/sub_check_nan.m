%%% Forward Euler checks
function bio = sub_check_nan(bio)
    %ID = (bio < 0);
    ID = (bio < 0 | isnan(bio));
    bio(ID) = eps();

end
