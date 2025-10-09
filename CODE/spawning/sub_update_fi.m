%%% Update biomass
function bio_out = sub_update_fi(bio_in,rec,nu,rep,gamma,die,nmort,fmort)
    % all inputs except rec & die are in g g-1 d-1; rec & die are g d-1
    % rec = rec from smaller size class = TOTAL biomass gained from recruitment
    % nu = energy avail for growth or repro
    % rep = energy lost to egg production
    % gamma = energy lost to maturation to larger size class
    % nmort = natural mortality rate
    % fmort = fishing rate
    % die = biomass lost to predation
    db = rec + ((nu - rep - gamma - nmort - fmort) .* bio_in) - die;
    bio_out = bio_in + db;
end
