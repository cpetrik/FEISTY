%%% Forward Euler checks
function bio = sub_check(bio)
    ID = (bio < 0);
    bio(ID) = eps();

end
