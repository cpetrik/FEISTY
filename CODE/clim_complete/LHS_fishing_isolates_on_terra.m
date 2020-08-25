% LHS fishing sets 
% Isolate runs with 2 types set at 0.2 or 0.3 and 3rd varied

load('/Volumes/FEISTY/NC/Matlab_new_size/bio_rates/LHS_diffF.mat');

%%
fid2 = find(X(:,1)==0.2);
fid3 = find(X(:,1)>0.29 & X(:,1)<0.31);

pid2 = find(X(:,2)==0.2);
pid3 = find(X(:,2)>0.29 & X(:,2)<0.31);

did2 = find(X(:,3)==0.2);
did3 = find(X(:,3)>0.29 & X(:,3)<0.31);

vP2 = intersect(fid2,did2);
vP3 = intersect(fid3,did3);

vD2 = intersect(fid2,pid2);
vD3 = intersect(fid3,pid3);

vF2 = intersect(pid2,did2);
vF3 = intersect(pid3,did3);

%%
X_new = X([vP2,vP3,vD2,vD3,vF2,vF3],:);

save('/Volumes/FEISTY/NC/Matlab_new_size/bio_rates/LHS_diffF3.mat','X_new');
