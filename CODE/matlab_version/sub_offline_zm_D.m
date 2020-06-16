%%% Offline coupling
function [out_1, zf] = sub_offline_zm_D(enc_1,bio_1,dZ)
    % ADD FLAG FOR COUNTING HOW MANY TIMES THIS HAPPENS
    % offline switch
    con_1 = enc_1 .* bio_1;
    
    out_1 = enc_1;
    
    id=((con_1) > dZ);
    
    frac1(id,1) = con_1(id,1) ./ (con_1(id,1));
    
    out_1(id,1) = (frac1(id,1) .* dZ(id,1)) ./ bio_1(id,1);
    
    zf = (out_1.*bio_1) ./ dZ;
end
