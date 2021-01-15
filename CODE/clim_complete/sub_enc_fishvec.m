%%%  Encounter rates
function enc = sub_enc_fishvec(param,A,fish)
    % bio: prey biomass
    % A: predator search rate
    % tprey: time spent in area with that prey item
    % theta: preference for prey item
    
    % Create tprey from td (tpel)
    tprey = repmat(fish.td,1,1,length(param.wc));
    %benthic prey are 1-td
    tprey(:,:,[3,10,11]) = 1 - tprey(:,:,[3,10,11]);
    
    enc = zeros(size(tprey));
    % Loop over each fn type
    for j = 1:length(param.ixFish)
        i = param.ixFish(j);
        enc(:,i,:) = fish.bio .* A(:,i) .* squeeze(tprey(:,i,:)) .* param.theta(i,:);
    end
    
    
end
