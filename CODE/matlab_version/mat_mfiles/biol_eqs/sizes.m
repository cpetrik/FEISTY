
%Zoop
L_zm = 10^((log10(0.2)+log10(2))/2);
L_zl = 10^((log10(2)+log10(20))/2);
M_zm = 9.0 * exp(1.953 + (2.399*log(L_zm)))*1.0e-6;
M_zl = 9.0 * exp(1.953 + (2.399*log(L_zl)))*1.0e-6;

%Fish
L_s = 10^((log10(2)+log10(20))/2); % small
L_m = 10^((log10(20)+log10(200))/2); % medium
L_l = 10^((log10(200)+log10(2000))/2); % large

%%! Mass from length using Andersen & Beyer 2013
% Convert from mm to cm and use their const coeff = 0.01g/cm3
M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

Ls=[L_s; 22/2; 10];
Lm=[L_m; 220/2; 200];
Ll=[L_l; 2200/2; 1000];
Ms = 0.01 .* (0.1.*Ls).^3;
Mm = 0.01 .* (0.1.*Lm).^3;
Ml = 0.01 .* (0.1.*Ll).^3;

%FishBase values
L_mf = 10*[45.6;52.5];
M_mf = [995;14815];
L_ld = 10*[135.1;82.3;58.2;106.5;83.1;92.1;120];
M_ld = [43157;18614;7808;14700;16544;19393;37490];
L_lp = 10*[112.1;151.7;221.4;122.2];
M_lp = [99640;9100;202690;48650];

mean(L_mf)
mean(M_mf)
mean(L_ld)
mean(M_ld)
mean(L_lp)
mean(M_lp)

%Work in opp direction from weight
W = [0.01;900;10e3];
L = (W./0.01) .^(1/3) ./0.1;
