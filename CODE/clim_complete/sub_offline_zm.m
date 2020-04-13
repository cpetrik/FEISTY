%%% Offline coupling
function [out_1, out_2, out_3, out_4, out_5, zf] = sub_offline_zm(enc_1,enc_2,enc_3,enc_4,enc_5,bio_1,bio_2,bio_3,bio_4,bio_5,dZ)
    % ADD FLAG FOR COUNTING HOW MANY TIMES THIS HAPPENS
    % offline switch
    con_1 = enc_1 .* bio_1;
    con_2 = enc_2 .* bio_2;
    con_3 = enc_3 .* bio_3;
    con_4 = enc_4 .* bio_4;
    con_5 = enc_5 .* bio_5;

    out_1 = enc_1;
    out_2 = enc_2;
    out_3 = enc_3;
    out_4 = enc_4;
    out_5 = enc_5;

    id=((con_1 + con_2 + con_3 + con_4 + con_5) > dZ);
    
    frac1(id,1) = con_1(id,1) ./ (con_1(id,1) + con_2(id,1) + con_3(id,1) + con_4(id,1) + con_5(id,1));
    frac2(id,1) = con_2(id,1) ./ (con_1(id,1) + con_2(id,1) + con_3(id,1) + con_4(id,1) + con_5(id,1));
    frac3(id,1) = con_3(id,1) ./ (con_1(id,1) + con_2(id,1) + con_3(id,1) + con_4(id,1) + con_5(id,1));
    frac4(id,1) = con_4(id,1) ./ (con_1(id,1) + con_2(id,1) + con_3(id,1) + con_4(id,1) + con_5(id,1));
    frac5(id,1) = con_5(id,1) ./ (con_1(id,1) + con_2(id,1) + con_3(id,1) + con_4(id,1) + con_5(id,1));
    
    out_1(id,1) = (frac1(id,1) .* dZ(id,1)) ./ bio_1(id,1);
    out_2(id,1) = (frac2(id,1) .* dZ(id,1)) ./ bio_2(id,1);
    out_3(id,1) = (frac3(id,1) .* dZ(id,1)) ./ bio_3(id,1);
    out_4(id,1) = (frac4(id,1) .* dZ(id,1)) ./ bio_4(id,1);
    out_5(id,1) = (frac5(id,1) .* dZ(id,1)) ./ bio_5(id,1);

    zf = (out_1.*bio_1 + out_2.*bio_2 + out_3.*bio_3 + out_4.*bio_4 + out_5.*bio_5) ./ dZ;
end
