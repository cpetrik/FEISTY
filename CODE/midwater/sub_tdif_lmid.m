%%% Fraction of time spent in pelagic (for demersal)
function tdif = sub_tdif_lmid(Z,param,mf,mm,mp,ml,md,bi)
    % bio1, bio2: pelagic prey
    % bio3, bio4: demersal prey

    % Varies with seafloor depth
    tdif = zeros(size(Z));
    cid = (Z < param.PI_be_cutoff);                           %continental shelf
    sid = (Z >= param.PI_be_cutoff & Z < param.MI_be_cutoff); %slope
    bid = (Z >= param.MI_be_cutoff);                          %basin
    % use preferences in calculation
    biop = zeros(size(Z));
    biom = zeros(size(Z));
    biod = zeros(size(Z));
    biop(cid) = param.LL_phi_MF(1) * mf(cid) + param.LL_phi_MP(1) * mp(cid);
    biom(cid) = param.LL_phi_MM(1) * mm(cid) + param.LL_phi_ML(1) * ml(cid);
    biod(cid) = param.LL_phi_MD(1) * md(cid) + param.LL_phi_BE(1) * bi(cid);
    biop(sid) = param.LL_phi_MF(2) * mf(sid) + param.LL_phi_MP(2) * mp(sid);
    biom(sid) = param.LL_phi_MM(2) * mm(sid) + param.LL_phi_ML(2) * ml(sid);
    biod(sid) = param.LL_phi_MD(2) * md(sid) + param.LL_phi_BE(2) * bi(sid);
    biop(bid) = param.LL_phi_MF(3) * mf(bid) + param.LL_phi_MP(3) * mp(bid);
    biom(bid) = param.LL_phi_MM(3) * mm(bid) + param.LL_phi_ML(3) * ml(bid);
    biod(bid) = param.LL_phi_MD(3) * md(bid) + param.LL_phi_BE(3) * bi(bid);


    tdiff(:,1) = biop ./ (biop + biom + biod);
    tdiff(:,2) = biom ./ (biop + biom + biod);
    tdiff(:,3) = biod ./ (biop + biom + biod);

    % Make sure it sums to one
    tdiff(:,1) = tdiff(:,1) ./ (tdiff(:,1) + tdiff(:,2) + tdiff(:,3));
    tdiff(:,2) = tdiff(:,2) ./ (tdiff(:,1) + tdiff(:,2) + tdiff(:,3));
    tdiff(:,3) = tdiff(:,3) ./ (tdiff(:,1) + tdiff(:,2) + tdiff(:,3));
    
end
