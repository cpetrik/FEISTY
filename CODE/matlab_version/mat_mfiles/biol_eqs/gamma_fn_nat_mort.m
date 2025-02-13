kap =1;
Z = 0.001/0.5;

fnu=mean(nSF(:));
pnu=mean(nSP(:)); 
dnu=mean(nSD(:));
nu = (fnu + pnu + dnu)/3;
D = [0.1:0.05:0.9]/365;

gg = ((kap.*nu) - D) ./ (1 - (Z.^(1 - (D ./ (kap.*nu)))));
gamma = min(gg,nu);