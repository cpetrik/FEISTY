
clearvars

% double click on input parameters and load 
% ben_prod, depthWOA, photic, RmaxL, mz_prod as a column

ben_prod = ben_prod(2:end);
depthWOA = depthWOA(2:end);
photic = photic(2:end);
RmaxL = RmaxL(2:end);
RmaxS = mz_prod(2:end); 

% remove data mz_prod

save globaldat_input.mat

clearvars

% load all temperature data as numeric
temp = inputparameters;
temp_grad = [0:5:100 125:25:500 550:50:2000 2100:100:5500];
temp_prof = temp(2:end,:);

% remove all except temp_grad and temp_prof
save globaldat_temp.mat