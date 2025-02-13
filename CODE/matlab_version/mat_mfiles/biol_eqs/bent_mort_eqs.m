% Test benthos carrying capacity

clear all
close all

load('/Volumes/GFDL/POEM_JLD/esm2m_hist/Data_ESM2Mhist_2005.mat');

ids = [40319,42639,41782,36334,38309,42744,30051,41284,38003,25248,33069];
ID = ids(1);

%%
BE=0.05;
bent=NaN*ones(366,1);
bent2=NaN*ones(366,1);
bent3=NaN*ones(366,1);
bent4=NaN*ones(366,1);
CC = 0.5;
eaten = rand(365,1)*5e-2;

%
bent(1)=1e-5;
bent2(1)=1e-5;
bent3(1)=1e-5;
bent4(1)=1e-5;
for t=1:365
    det = COBALT.det(ID,t);
    r1 = BE*det;%/bent(t);
    r2 = BE*det/bent2(t);
    r3 = BE*det;
    %eaten = 4e-3; %rand; %con.*bio;
    bent(t+1) = max(eps, bent(t) + r1 * (CC - bent(t)) - eaten(t)); %chemostat
    bent2(t+1) = max(eps, bent2(t) + r2 * bent2(t) * (1 - bent2(t)/CC) - eaten(t)); %logistic
    bent3(t+1) = max(eps, bent3(t) + r3 - eaten(t)); %old
    bent4(t+1) = max(eps,bent4(t) + BE*det * (1 - bent4(t)/CC) - eaten(t)); 
    
end

test = BE*COBALT.det(ID,:);
%%
figure
plot(bent,'b','LineWidth',2); hold on;
plot(bent2,'--r','LineWidth',2); hold on;
plot(bent3,'-.g','LineWidth',2); hold on;
plot(bent4,':m','LineWidth',2); hold on;
plot(test,'k','LineWidth',2); hold on;
ylim([0 1])