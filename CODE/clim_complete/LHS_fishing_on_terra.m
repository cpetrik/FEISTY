%LHS fishing sets remaining

load('/Volumes/FEISTY/NC/Matlab_new_size/bio_rates/LHS_diffF.mat');

%%
nw = floor(length(X_rem)/10); 

wsets = 1:nw:length(X_rem);
wcomp = wsets+49;
wend = X_rem(wcomp(1:10),:);

%%
mis2=[];
for i=1:10
    mis2 = [mis2 wcomp(i):wsets(i+1)];
end

X_rem2 = X_rem(mis2,:);

save('/Volumes/FEISTY/NC/Matlab_new_size/bio_rates/LHS_diffF2.mat');