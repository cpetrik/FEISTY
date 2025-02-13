%%% Offline coupling
function [bf] = sub_offline_bent(con,bio,B)
    tcon = sum(con .* bio,2);
    bf = tcon ./ B; 
end
