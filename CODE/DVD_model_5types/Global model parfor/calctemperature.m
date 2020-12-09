
% calculate temperature dependencies integrating the vertical distribution

function[scTemp, scTempm] = calctemperature(param,gridi,temp_grad,temp_prof,depthDay,depthNight,bottom)

tempgrid = temp_prof(gridi,:); % get temperature profile of the grid
tempdat = zeros(bottom+1,1); % get vector for each meter

clo = abs(temp_grad - bottom); % find index that corresponds to seafloor
minDist = min(clo);
idxn = find(clo == minDist);

for i = 1:(idxn-1)
    temp = linspace(tempgrid(i),tempgrid(i+1),(temp_grad(i+1)-temp_grad(i)+1));
    tempdat((temp_grad(i)+1):temp_grad(i+1)) = temp(1:(end-1));
end
tempdat(end) = tempgrid(idxn);

dist = (depthDay + depthNight)/2;
TQ10 =  param.Q10.^((tempdat(1:bottom+1)-10)/10);
TQ10m =  param.Q10m.^((tempdat(1:bottom+1)-10)/10);

scTemp_step = dist .* TQ10;
scTemp = sum(scTemp_step,1);

scTemp_stepm = dist .* TQ10m;
scTempm = sum(scTemp_stepm,1);

