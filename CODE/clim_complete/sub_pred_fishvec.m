%%%  Predation rates
function pred = sub_pred_fishvec(param,fish)
    % bio: prey biomass
    
    pred = zeros(param.NX,length(param.wc));
    % Loop over each fn type
    for j = 1:length(param.ixFish)
        i = param.ixFish(j);
        %tot pred = sum(predator biom * predator con)
        pred(:,i) = sum(fish.bio .* squeeze(fish.con(:,:,i)),2);
    end
    
end
