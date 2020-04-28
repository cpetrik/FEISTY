%%% Offline coupling
function [bf] = sub_offline_be2(con,bio,Bs,Bm)
    tcon = sum(con .* bio,2);
    B = Bs + Bm;
    bf = tcon ./ B; 
end
