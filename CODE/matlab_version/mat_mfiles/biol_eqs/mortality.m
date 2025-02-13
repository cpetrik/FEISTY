% Mortality rates by size

clear all
close all

L_s = 10^((log10(2)+log10(20))/2);
L_m = 10^((log10(20)+log10(200))/2);
L_l = 10^((log10(200)+log10(2000))/2);

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;
s = [M_s;M_m;M_l];

temp=-2:30;

M = repmat(0.44,size(s));
AB = 0.35 .* 4.5 .* s.^(-0.25);
H = 0.84 .* s.^(-0.25);

%%
figure
subplot(2,2,1)
plot(log10(s),M,'.k','MarkerSize',25); hold on;
plot(log10(s),AB,'.b','MarkerSize',25); hold on;
plot(log10(s),H,'.r','MarkerSize',25); hold on;
%ylim([0 2])
legend('M','AB','H')