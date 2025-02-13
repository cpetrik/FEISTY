% Kiorboe & Hirst eqs
% Play with clearance rate for MF

clear all
close all

L_s = 10^((log10(2)+log10(20))/2);
L_m = 10^((log10(20)+log10(200))/2);
L_l = 10^((log10(200)+log10(2000))/2);

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

%% Kiorboe & Hirst orig
% T=15C
temp=15;

%Max clearance rate/search volume
FmaxLarv = exp(0.063*(temp-15.0)) .* 10^(8.1) .* M_s^(0.057) .* (24e-3/9);
FmaxS = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_s^(-0.24) .* (24e-3/9);
FmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_m^(-0.24) .* (24e-3/9);
FmaxL = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_l^(-0.24) .* (24e-3/9);

m=[M_s; M_m; M_l];
fmax=[FmaxS;FmaxM;FmaxL];

figure(1)
%subplot(2,2,1)
plot(log10(m),log10(fmax),'.b','MarkerSize',25); hold on;
plot(log10(M_s),log10(FmaxLarv),'.m','MarkerSize',25); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H','K&H larv')
legend('location','southwest')
title('Specific clearance rate')
print('-dpng','Clearance_mass_orig.png')


%% vs temp
temp=-2:30;

%Max clearance rate/search volume
FmaxLarv = exp(0.063*(temp-15.0)) .* 10^(8.1) .* M_s^(0.057) .* (24e-3/9);
FmaxS = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_s^(-0.24) .* (24e-3/9);
FmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_m^(-0.24) .* (24e-3/9);
FmaxL = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_l^(-0.24) .* (24e-3/9);

figure(2)
%subplot(2,2,1)
plot(temp,log10(FmaxS),'r','LineWidth',2); hold on;
plot(temp,log10(FmaxM),'k','LineWidth',2); hold on;
plot(temp,log10(FmaxL),'b','LineWidth',2); hold on;
plot(temp,log10(FmaxLarv),'m','LineWidth',2); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('temp (C)')
legend('K&H S','K&H M','K&H L','K&H larv')
legend('location','west')
title('Specific clearance rate')
print('-dpng','Clearance_temp_orig.png')


%% Kiorboe & Hirst NEW Larvae
% T=15C
temp=15;

%Max clearance rate/search volume
FmaxLarv = exp(0.063*(temp-15.0)) .* 10^(3.58) .* M_s^(0.057) .* (24e-3/9);
FmaxS = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_s^(-0.24) .* (24e-3/9);
FmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_m^(-0.24) .* (24e-3/9);
FmaxL = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_l^(-0.24) .* (24e-3/9);

m=[M_s; M_m; M_l];
fmax=[FmaxS;FmaxM;FmaxL];

figure(3)
%subplot(2,2,1)
plot(log10(m),log10(fmax),'.b','MarkerSize',25); hold on;
plot(log10(M_s),log10(FmaxLarv),'.m','MarkerSize',25); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H','K&H larv')
legend('location','southwest')
title('Specific clearance rate')
print('-dpng','Clearance_mass_newLarv.png')


%% vs temp
temp=-2:30;

%Max clearance rate/search volume
FmaxLarv = exp(0.063*(temp-15.0)) .* 10^(3.58) .* M_s^(0.057) .* (24e-3/9);
FmaxS = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_s^(-0.24) .* (24e-3/9);
FmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_m^(-0.24) .* (24e-3/9);
FmaxL = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_l^(-0.24) .* (24e-3/9);

figure(4)
%subplot(2,2,1)
plot(temp,log10(FmaxS),'r','LineWidth',2); hold on;
plot(temp,log10(FmaxM),'k','LineWidth',2); hold on;
plot(temp,log10(FmaxL),'b','LineWidth',2); hold on;
plot(temp,log10(FmaxLarv),'m','LineWidth',2); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('temp (C)')
legend('K&H S','K&H M','K&H L','K&H larv')
legend('location','west')
title('Specific clearance rate')
print('-dpng','Clearance_comp_temp_newLarv.png')


%% Kiorboe & Hirst NEW All
% T=15C
temp=15;

%Max clearance rate/search volume
FmaxLarv = exp(0.063*(temp-15.0)) .* 10^(3.68) .* M_s^(0.1963) .* (24e-3/9);
FmaxS = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_s^(-0.24) .* (24e-3/9);
FmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_m^(-0.24) .* (24e-3/9);
FmaxL = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_l^(-0.24) .* (24e-3/9);
FmaxMF1 = exp(0.063*(temp-15.0)) .* 10^(3.44) .* M_m^(-0.29) .* (24e-3/9);
FmaxMP1 = exp(0.063*(temp-15.0)) .* 10^(2.87) .* M_m^(-0.29) .* (24e-3/9);
FmaxMF2 = exp(0.063*(temp-15.0)) .* 10^(3.60) .* M_m^(-0.36) .* (24e-3/9);
FmaxMP2 = exp(0.063*(temp-15.0)) .* 10^(3.13) .* M_m^(0.06) .* (24e-3/9);

m=[M_s; M_m; M_l];
fmax=[FmaxS;FmaxM;FmaxL];

figure(5)
%subplot(2,2,1)
plot(log10(m),log10(fmax),'.k','MarkerSize',25); hold on;
plot(log10(M_s),log10(FmaxLarv),'.m','MarkerSize',25); hold on;
plot(log10(M_m),log10(FmaxMF1),'.b','MarkerSize',25); hold on;
plot(log10(M_m),log10(FmaxMP1),'.r','MarkerSize',25); hold on;
plot(log10(M_m),log10(FmaxMF2),'*b','MarkerSize',5); hold on;
plot(log10(M_m),log10(FmaxMP2),'*r','MarkerSize',5); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H','K&H larv','Clupe1','Other1','Clupe2','Other2')
legend('location','southwest')
title('Specific clearance rate')
print('-dpng','Clearance_mass_new.png')


%% vs mass
temp=15;
lW=-4:5;
W=10.^(lW);

%Max clearance rate/search volume
FmaxLarv = exp(0.063*(temp-15.0)) .* 10^(3.68) .* W.^(0.1963) .* (24e-3/9);
FmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* W.^(-0.24) .* (24e-3/9);
FmaxMF1 = exp(0.063*(temp-15.0)) .* 10^(3.44) .* W.^(-0.29) .* (24e-3/9);
FmaxMP1 = exp(0.063*(temp-15.0)) .* 10^(2.87) .* W.^(-0.29) .* (24e-3/9);
FmaxMF2 = exp(0.063*(temp-15.0)) .* 10^(3.60) .* W.^(-0.36) .* (24e-3/9);
FmaxMP2 = exp(0.063*(temp-15.0)) .* 10^(3.13) .* W.^(0.06) .* (24e-3/9);

figure(6)
%subplot(2,2,1)
plot(lW,log10(FmaxM),'k','LineWidth',2); hold on;
plot(lW,log10(FmaxMF1),'c','LineWidth',2); hold on;
plot(lW,log10(FmaxMP1),'g','LineWidth',2); hold on;
plot(lW,log10(FmaxMF2),'c--','LineWidth',2); hold on;
plot(lW,log10(FmaxMP2),'g--','LineWidth',2); hold on;
plot(lW,log10(FmaxLarv),'m--','LineWidth',2); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H','Clupe1','Other1','Clupe2','Other2','K&H larv')
legend('location','southwest')
title('Specific clearance rate')
print('-dpng','Clearance_mass_fn_new.png')


%% vs temp
temp=-2:30;

%Max clearance rate/search volume
FmaxLarv = exp(0.063*(temp-15.0)) .* 10^(3.68) .* M_s^(0.1963) .* (24e-3/9);
FmaxS = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_s^(-0.24) .* (24e-3/9);
FmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_m^(-0.24) .* (24e-3/9);
FmaxL = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_l^(-0.24) .* (24e-3/9);
FmaxMF1 = exp(0.063*(temp-15.0)) .* 10^(3.44) .* M_m^(-0.29) .* (24e-3/9);
FmaxMP1 = exp(0.063*(temp-15.0)) .* 10^(2.87) .* M_m^(-0.29) .* (24e-3/9);
FmaxMF2 = exp(0.063*(temp-15.0)) .* 10^(3.60) .* M_m^(-0.36) .* (24e-3/9);
FmaxMP2 = exp(0.063*(temp-15.0)) .* 10^(3.13) .* M_m^(0.06) .* (24e-3/9);

figure(7)
%subplot(2,2,1)
plot(temp,log10(FmaxS),'r','LineWidth',2); hold on;
plot(temp,log10(FmaxM),'k','LineWidth',2); hold on;
plot(temp,log10(FmaxMF1),'c','LineWidth',2); hold on;
plot(temp,log10(FmaxMP1),'g','LineWidth',2); hold on;
plot(temp,log10(FmaxMF2),'c--','LineWidth',2); hold on;
plot(temp,log10(FmaxMP2),'g--','LineWidth',2); hold on;
plot(temp,log10(FmaxL),'b','LineWidth',2); hold on;
plot(temp,log10(FmaxLarv),'m--','LineWidth',2); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('temp (C)')
legend('K&H S','K&H M','Clupe1','Other1','Clupe2','Other2','K&H L','K&H larv')
legend('location','northwest')
title('Specific clearance rate')
print('-dpng','Clearance_comp_temp_new.png')

%% Figures of just MF values used in simulations

temp=15;
%Max clearance rate/search volume
sFmaxS = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_s^(-0.24) .* (24e-3/9);
sFmaxSF1 = exp(0.063*(temp-15.0)) .* 10^(3.44) .* M_s^(-0.29) .* (24e-3/9);
sFmaxSF2 = exp(0.063*(temp-15.0)) .* 10^(3.60) .* M_s^(-0.36) .* (24e-3/9);
sFmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_m^(-0.24) .* (24e-3/9);
sFmaxMF1 = exp(0.063*(temp-15.0)) .* 10^(3.44) .* M_m^(-0.29) .* (24e-3/9);
sFmaxMF2 = exp(0.063*(temp-15.0)) .* 10^(3.60) .* M_m^(-0.36) .* (24e-3/9);

mFmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* W.^(-0.24) .* (24e-3/9);
mFmaxMF1 = exp(0.063*(temp-15.0)) .* 10^(3.44) .* W.^(-0.29) .* (24e-3/9);
mFmaxMF2 = exp(0.063*(temp-15.0)) .* 10^(3.60) .* W.^(-0.36) .* (24e-3/9);

temp=-2:30;
tFmaxM = exp(0.063*(temp-15.0)) .* 1.74e3 .* M_m^(-0.24) .* (24e-3/9);
tFmaxMF1 = exp(0.063*(temp-15.0)) .* 10^(3.44) .* M_m^(-0.29) .* (24e-3/9);
tFmaxMF2 = exp(0.063*(temp-15.0)) .* 10^(3.60) .* M_m^(-0.36) .* (24e-3/9);

figure(8)
subplot(2,2,1)
plot(log10(M_s),log10(sFmaxS),'.k','MarkerSize',25); hold on;
plot(log10(M_s),log10(sFmaxSF1),'.b','MarkerSize',25); hold on;
plot(log10(M_s),log10(sFmaxSF2),'.r','MarkerSize',25); hold on;
plot(log10(M_m),log10(sFmaxM),'.k','MarkerSize',25); hold on;
plot(log10(M_m),log10(sFmaxMF1),'.b','MarkerSize',25); hold on;
plot(log10(M_m),log10(sFmaxMF2),'.r','MarkerSize',25); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H','F1','F2','K&H','F1','F2')
legend('location','northeast')
title('Specific clearance rate')
xlim([-4 5])

subplot(2,2,3)
plot(lW,log10(mFmaxM),'k','LineWidth',2); hold on;
plot(lW,log10(mFmaxMF1),'b','LineWidth',2); hold on;
plot(lW,log10(mFmaxMF2),'r','LineWidth',2); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H all','Clupe 1','Clupe 2')
legend('location','southwest')
title('Specific clearance rate')
xlim([-4 5])

subplot(2,2,2)
plot(temp,log10(tFmaxM),'k','LineWidth',2); hold on;
plot(temp,log10(tFmaxMF1),'b','LineWidth',2); hold on;
plot(temp,log10(tFmaxMF2),'r','LineWidth',2); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('temp (C)')
legend('K&H all','MF1','MF2')
legend('location','northwest')
title('Specific clearance rate')
xlim([-2 30])
print('-dpng','Clearance_new_MF.png')