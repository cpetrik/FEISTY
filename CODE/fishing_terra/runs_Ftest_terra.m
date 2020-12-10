load('/Volumes/FEISTY/NC/Matlab_new_size/bio_rates/LHS_diffF.mat','X');


terra_st = [1,908,1815,2721,3627,4533,5439,6345,7251,8157] + 199;
terra_en = terra_st + 49;
Compl = 1:200;
for i=1:length(terra_st)
    Compl = [Compl terra_st(i):terra_en(i)];
end

Compl = Compl';
all = 1:length(X);
mis = setdiff(all,Compl);

X_rem = X(mis,:);

save('/Volumes/FEISTY/NC/Matlab_new_size/bio_rates/LHS_diffF.mat','X_rem',...
    'Compl','mis','-append');