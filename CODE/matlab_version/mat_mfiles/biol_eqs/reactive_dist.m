% Reactive distance

%Most are from larval studies
L_zm = 10^((log10(0.2)+log10(2))/2);
L_s = 10^((log10(2)+log10(20))/2); 
L_m = 10^((log10(20)+log10(200))/2); 
L_l = 10^((log10(200)+log10(2000))/2); 
%Prey lengths: Zm, Zl/S, M, L
L=[L_zm;L_s;L_m];
%Predator lengths: S, M, L
BL=[L_s;L_m;L_l];

%% r^2 * exp(c*r) = (rho * E0 * exp(-Kz) * C * pi * B^2) / se
C=0.5;
E0=100;
K=1.14;
c=1.57;
rho=0.5;
z=100;

se = 1.41e-6 ./ (0.8 + 0.03.*BL).^3;
r1 = (rho .* E0 .* exp(-K.*z) .* C .* pi .* L.^2) ./ se;

figure(1)
plot(BL,r1)

%% r^2 * exp(c*r) = C * A * E * (eb/(ke+eb));
C=0.3;
A=pi.*(0.5.*L).^2;
ke=1;
eb=100;
Eh = 10.^ (4.88 ./ (1 + exp(-(BL - 10.98) ./ 1.34)));
Ec = 10.^ (5.04 ./ (1 + exp(-(BL - 5.33) ./ 0.64)));

rher = sqrt(C .* A .* Eh .* (eb./(ke+eb))); %herring
rcod = sqrt(C .* A .* Ec .* (eb./(ke+eb))); %cod

rp = exp(-0.163 + 0.717 .* log(BL)); %pause-travel
rc = exp(-0.73 + 0.931 .* log(BL)); %cruising

x=0:BL(end);
y=0:BL(end);
j=3*y;

%%
figure(2)
plot(BL,rher,'r','LineWidth',2); hold on;
plot(BL,rcod,'b','LineWidth',2); hold on;
plot(BL,rp,'c','LineWidth',2); hold on;
plot(BL,rc,'k','LineWidth',2); hold on;
plot(x,j,'m','LineWidth',2); hold on;
plot(x,y,'--k');

