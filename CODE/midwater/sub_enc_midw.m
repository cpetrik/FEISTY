%%%  Encounter rates
function enc = sub_enc_midw(H,param,Tp,Tb,wgt,prey,tdif,tprey,pref)
    %H: depth
    %Tp: epipel temp
    %Tm: mesopel temp
    %Tb: bottom temp
    % wgt: ind weight of size class
    % pred: pred biomass density,
    % prey: prey biomass density,
    % A: predator search rate,
    % tdif: frac time (1) epi, (2) meso, (3) benthic - for calculating mean temp
    % tprey: time spent in area with that prey item - for calculating enc
    % pref: preference for prey item

    temp = (Tp.*tdif(:,1)) + (Tm.*tdif(:,2)) + (Tb.*tdif(:,3));

    %! Clearance rate
    A = (exp(param.ke * (temp-10.0)) .* param.gam .* wgt^(-param.benc)) ./365.0;

    %! Encounter per predator, mult by biomass later

    %Prefs vary by depth
    cid = (H < param.PI_be_cutoff);                           %continental shelf
    sid = (H >= param.PI_be_cutoff & H < param.MI_be_cutoff); %slope
    bid = (H >= param.MI_be_cutoff);                          %basin

    %Overlap of predator and prey = tdif * tprey
    % frac = zeros(param.NX,1);
    % ID = (tprey>0);
    % frac(ID) = 1.0;
    frac = tprey;

    enc = zeros(param.NX,1);
    enc(cid,:) = prey(cid,:) .* A(cid,:) .* frac .* pref(1);
    enc(sid,:) = prey(sid,:) .* A(sid,:) .* frac .* pref(2);
    enc(bid,:) = prey(bid,:) .* A(bid,:) .* frac .* pref(3);
end
