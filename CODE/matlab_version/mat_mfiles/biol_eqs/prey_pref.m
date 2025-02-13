% Look at prey preference with James' old code

clear all
close all

L_zm = 10^((log10(0.2)+log10(2))/2);
L_zl = 10^((log10(2)+log10(20))/2);
L_s = 10^((log10(2)+log10(20))/2);
L_m = 10^((log10(20)+log10(200))/2);
L_l = 10^((log10(200)+log10(2000))/2);

M_zm = 9.0 * exp(1.953 + (2.399*log(L_zm))) *1.0e-6;
M_zl = 9.0 * exp(1.953 + (2.399*log(L_zl))) *1.0e-6;
M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

%%
% w1: body size of pred (g)
w1 = [M_s;M_m;M_l];
% w2: body size of prey (g)
w2 = [M_zm;M_zl;M_s;M_m;M_l];
beta = 2; % mean PPMR
sigma = 1; % sd PPMR

PIJ = zeros(length(w2),length(w1));
for j=1:length(w1)
    for i = 1:length(w2)
        x1 = log10(w1(j));
        x2 = log10(w2(i));
        if w1(j) > w2(i)
            pij = exp(-((x1-x2-beta).^2) / (2.*sigma.^2));
        else
            pij = 0;
        end
        PIJ(i,j) = pij;
    end
end


