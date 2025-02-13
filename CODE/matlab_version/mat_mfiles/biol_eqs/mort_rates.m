%Mortality rates

clear all
close all

L_s = 10^((log10(2)+log10(20))/2);
L_m = 10^((log10(20)+log10(200))/2);
L_l = 10^((log10(200)+log10(2000))/2);

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;
M=[M_s;M_m;M_l];

temp=-2:30;

ab(1,1) = 0.35 * 4.5 * M_s^(-0.25); %(includes predation, excludes fishing)
ab(2,1) = 0.35 * 4.5 * M_m^(-0.25);
ab(3,1) = 0.35 * 4.5 * M_l^(-0.25);

Nat = (0.44 / 365) * ones(3,1);

abT(1,:) = 0.35 * 4.5 * M_s.^(-0.25) .* exp(0.063.*(temp-15.0));
abT(2,:) = 0.35 * 4.5 * M_m.^(-0.25) .* exp(0.063.*(temp-15.0));
abT(3,:) = 0.35 * 4.5 * M_l.^(-0.25) .* exp(0.063.*(temp-15.0));
NatT = (0.44 / 365) .* exp(0.063.*(temp-15.0));

%%
figure(1)
plot(log(M),ab/365,'b.','MarkerSize',20); hold on;
plot(log(M),Nat,'--k','LineWidth',2)
xlabel('ln Mass')
ylabel('Daily mortality')
legend('AB','Megrey')

figure(2)
plot(temp,abT/365,'LineWidth',2); hold on;
plot(temp,NatT,'k','LineWidth',2)
xlabel('Temp')
ylabel('Daily mortality')
legend('ABs','ABm','ABl','Megrey')
