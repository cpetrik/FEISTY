clear all
close all

%!Individual Mass (g)
M_s = 10^((log10(0.001)+log10(0.5))/2);
M_m = 10^((log10(0.5)+log10(250))/2);
M_l = 10^((log10(250)+log10(125000))/2);
m=[M_s; M_m; M_l];

kt=0.0855;
amet = 4;
bpow = 0.175;
T = 10;

% Metabolism
met15 = (exp(kt*(T-10.0)) .* amet .* m.^(-0.15)) ./365.0;
met17 = (exp(kt*(T-10.0)) .* amet .* m.^(-0.175)) ./365.0;
met20 = (exp(kt*(T-10.0)) .* amet .* m.^(-0.2)) ./365.0;

%Enc
A15 = (exp(0.063*(T-10.0)) .* 70 .* m.^(-0.15)) ./365.0;
A17 = (exp(0.063*(T-10.0)) .* 70 .* m.^(-0.175)) ./365.0;
A20 = (exp(0.063*(T-10.0)) .* 70 .* m.^(-0.20)) ./365.0;

%%
figure(1)
plot(log(m),log(met15),'r.-'); hold on;
plot(log(m),log(met17),'b.-'); hold on;
plot(log(m),log(met20),'k.-'); hold on;
legend('15','17','20')

figure(2)
plot(log(m),log(A15),'r.-'); hold on;
plot(log(m),log(A17),'b.-'); hold on;
plot(log(m),log(A20),'k.-'); hold on;
legend('15','17','20')
