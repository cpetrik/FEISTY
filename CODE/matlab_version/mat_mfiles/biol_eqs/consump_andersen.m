% Ken Andersen consumption 

h=60.0 / 365.0;
flev1=3500;
flev2=4500;
flev3=5500;
q=1;

L_s = 10^((log10(2)+log10(20))/2); 
L_m = 10^((log10(20)+log10(200))/2); 
L_l = 10^((log10(200)+log10(2000))/2); 

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

T = -2:0.2:37;

R=0:10;

%Andersen
Acmax_s = exp(0.063*(T-15.0)) * h * M_s^(3/4);
Acmax_m = exp(0.063*(T-15.0)) * h * M_m^(3/4);
Acmax_l = exp(0.063*(T-15.0)) * h * M_l^(3/4);

beta_s1 = exp(0.063*T-15.0) * flev1 * M_s^(q);
beta_m1 = exp(0.063*T-15.0) * flev1 * M_m^(q);
beta_l1 = exp(0.063*T-15.0) * flev1 * M_l^(q);
beta_s2 = exp(0.063*T-15.0) * flev2 * M_s^(q);
beta_m2 = exp(0.063*T-15.0) * flev2 * M_m^(q);
beta_l2 = exp(0.063*T-15.0) * flev2 * M_l^(q);
beta_s3 = exp(0.063*T-15.0) * flev3 * M_s^(q);
beta_m3 = exp(0.063*T-15.0) * flev3 * M_m^(q);
beta_l3 = exp(0.063*T-15.0) * flev3 * M_l^(q);

%%
figure(1)
subplot(3,1,1)
plot(T,Acmax_s,'k','LineWidth',2); hold on;
title('Small/Larvae 6mm 2.5mg')
subplot(3,1,2)
plot(T,Acmax_m,'k','LineWidth',2); hold on;
ylabel('Cmax (g d^-^1)')
title('Medium/Juveniles 63mm 2.5g')
subplot(3,1,3)
plot(T,Acmax_l,'k','LineWidth',2); hold on;
title('Large/Adults 630mm 2.5kg')
xlabel('Temp')
print('-dpng','Cmax_vs_temp.png')

figure(2)
subplot(3,1,1)
plot(T,beta_s1,'b','LineWidth',2); hold on;
plot(T,beta_s2,'r','LineWidth',2); hold on;
plot(T,beta_s3,'k','LineWidth',2); hold on;
title('Small/Larvae 6mm 2.5mg')
legend('g=3500','g=4500','g=5500')
legend('location','northwest')
subplot(3,1,2)
plot(T,beta_m1,'b','LineWidth',2); hold on;
plot(T,beta_m2,'r','LineWidth',2); hold on;
plot(T,beta_m3,'k','LineWidth',2); hold on;
title('Medium/Juveniles 63mm 2.5g')
ylabel('Beta (d^-^1)')
subplot(3,1,3)
plot(T,beta_l1,'b','LineWidth',2); hold on;
plot(T,beta_l2,'r','LineWidth',2); hold on;
plot(T,beta_l3,'k','LineWidth',2); hold on;
title('Large/Adults 630mm 2.5kg')
xlabel('Temp')
print('-dpng','Beta_vs_temp.png')

%%
%at T-10C
Acon_s1 = Acmax_s(61) .* (beta_s1(61) .* R) ./ (Acmax_s(61) + beta_s1(61).*R);
Acon_m1 = Acmax_m(61) .* (beta_m1(61) .* R) ./ (Acmax_m(61) + beta_s1(61).*R);
Acon_l1 = Acmax_l(61) .* (beta_l1(61) .* R) ./ (Acmax_l(61) + beta_s1(61).*R);
Acon_s2 = Acmax_s(61) .* (beta_s2(61) .* R) ./ (Acmax_s(61) + beta_s2(61).*R);
Acon_m2 = Acmax_m(61) .* (beta_m2(61) .* R) ./ (Acmax_m(61) + beta_s2(61).*R);
Acon_l2 = Acmax_l(61) .* (beta_l2(61) .* R) ./ (Acmax_l(61) + beta_s2(61).*R);
Acon_s3 = Acmax_s(61) .* (beta_s3(61) .* R) ./ (Acmax_s(61) + beta_s3(61).*R);
Acon_m3 = Acmax_m(61) .* (beta_m3(61) .* R) ./ (Acmax_m(61) + beta_s3(61).*R);
Acon_l3 = Acmax_l(61) .* (beta_l3(61) .* R) ./ (Acmax_l(61) + beta_s3(61).*R);


Acon_s4 = Acmax_s(61) .* (1 .* R) ./ (Acmax_s(61) + 1.*R);
Acon_m4 = Acmax_m(61) .* (1 .* R) ./ (Acmax_m(61) + 1.*R);
Acon_l4 = Acmax_l(61) .* (1 .* R) ./ (Acmax_l(61) + 1.*R);

%%
figure(3)
subplot(3,1,1)
plot(R,Acon_s1,'b','LineWidth',2); hold on;
plot(R,Acon_s2,'r','LineWidth',2); hold on;
plot(R,Acon_s3,'k','LineWidth',2); hold on;
%plot(R,Acon_s4,'m','LineWidth',2); hold on;
subplot(3,1,2)
plot(R,Acon_m1,'b','LineWidth',2); hold on;
plot(R,Acon_m2,'r','LineWidth',2); hold on;
plot(R,Acon_m3,'k','LineWidth',2); hold on;
%plot(R,Acon_m4,'m','LineWidth',2); hold on;
subplot(3,1,3)
plot(R,Acon_l1,'b','LineWidth',2); hold on;
plot(R,Acon_l2,'r','LineWidth',2); hold on;
plot(R,Acon_l3,'k','LineWidth',2); hold on;
%plot(R,Acon_l4,'m','LineWidth',2); hold on;

%% Feeding level
flev1=3500;
flev2=4500;
flev3=5500;

T=5;
Acmax_s = exp(0.063*(T-15.0)) * h * M_s^(3/4);
Acmax_m = exp(0.063*(T-15.0)) * h * M_m^(3/4);
Acmax_l = exp(0.063*(T-15.0)) * h * M_l^(3/4);
beta_s1 = exp(0.063*T-15.0) * flev1 * M_s^(q);
beta_m1 = exp(0.063*T-15.0) * flev1 * M_m^(q);
beta_l1 = exp(0.063*T-15.0) * flev1 * M_l^(q);
beta_s2 = exp(0.063*T-15.0) * flev2 * M_s^(q);
beta_m2 = exp(0.063*T-15.0) * flev2 * M_m^(q);
beta_l2 = exp(0.063*T-15.0) * flev2 * M_l^(q);
beta_s3 = exp(0.063*T-15.0) * flev3 * M_s^(q);
beta_m3 = exp(0.063*T-15.0) * flev3 * M_m^(q);
beta_l3 = exp(0.063*T-15.0) * flev3 * M_l^(q);

f_s1 = beta_s1.*R ./ (Acmax_s + beta_s1.*R);
f_s2 = beta_s2.*R ./ (Acmax_s + beta_s2.*R);
f_s3 = beta_s3.*R ./ (Acmax_s + beta_s3.*R);
f_m1 = beta_m1.*R ./ (Acmax_m + beta_m1.*R);
f_m2 = beta_m2.*R ./ (Acmax_m + beta_m2.*R);
f_m3 = beta_m3.*R ./ (Acmax_m + beta_m3.*R);
f_l1 = beta_l1.*R ./ (Acmax_l + beta_l1.*R);
f_l2 = beta_l2.*R ./ (Acmax_l + beta_l2.*R);
f_l3 = beta_l3.*R ./ (Acmax_l + beta_l3.*R);

figure(4)
subplot(3,1,1)
plot(R,f_s1,'b','LineWidth',2); hold on;
plot(R,f_s2,'r','LineWidth',2); hold on;
plot(R,f_s3,'k','LineWidth',2); hold on;
title('T=5C')
legend('g=3500','g=4500','g=5500')
subplot(3,1,2)
plot(R,f_m1,'b','LineWidth',2); hold on;
plot(R,f_m2,'r','LineWidth',2); hold on;
plot(R,f_m3,'k','LineWidth',2); hold on;
ylabel('Feeding level')
subplot(3,1,3)
plot(R,f_l1,'b','LineWidth',2); hold on;
plot(R,f_l2,'r','LineWidth',2); hold on;
plot(R,f_l3,'k','LineWidth',2); hold on;
xlabel('prey density')


T=15;
Acmax_s = exp(0.063*(T-15.0)) * h * M_s^(3/4);
Acmax_m = exp(0.063*(T-15.0)) * h * M_m^(3/4);
Acmax_l = exp(0.063*(T-15.0)) * h * M_l^(3/4);
beta_s1 = exp(0.063*T-15.0) * flev1 * M_s^(q);
beta_m1 = exp(0.063*T-15.0) * flev1 * M_m^(q);
beta_l1 = exp(0.063*T-15.0) * flev1 * M_l^(q);
beta_s2 = exp(0.063*T-15.0) * flev2 * M_s^(q);
beta_m2 = exp(0.063*T-15.0) * flev2 * M_m^(q);
beta_l2 = exp(0.063*T-15.0) * flev2 * M_l^(q);
beta_s3 = exp(0.063*T-15.0) * flev3 * M_s^(q);
beta_m3 = exp(0.063*T-15.0) * flev3 * M_m^(q);
beta_l3 = exp(0.063*T-15.0) * flev3 * M_l^(q);

f_s1 = beta_s1.*R ./ (Acmax_s + beta_s1.*R);
f_s2 = beta_s2.*R ./ (Acmax_s + beta_s2.*R);
f_s3 = beta_s3.*R ./ (Acmax_s + beta_s3.*R);
f_m1 = beta_m1.*R ./ (Acmax_m + beta_m1.*R);
f_m2 = beta_m2.*R ./ (Acmax_m + beta_m2.*R);
f_m3 = beta_m3.*R ./ (Acmax_m + beta_m3.*R);
f_l1 = beta_l1.*R ./ (Acmax_l + beta_l1.*R);
f_l2 = beta_l2.*R ./ (Acmax_l + beta_l2.*R);
f_l3 = beta_l3.*R ./ (Acmax_l + beta_l3.*R);

figure(5)
subplot(3,1,1)
plot(R,f_s1,'b','LineWidth',2); hold on;
plot(R,f_s2,'r','LineWidth',2); hold on;
plot(R,f_s3,'k','LineWidth',2); hold on;
title('Small/Larvae 6mm 2.5mg')
legend('g=3500','g=4500','g=5500')
legend('location','northwest')
subplot(3,1,2)
plot(R,f_m1,'b','LineWidth',2); hold on;
plot(R,f_m2,'r','LineWidth',2); hold on;
plot(R,f_m3,'k','LineWidth',2); hold on;
title('Medium/Juveniles 63mm 2.5g')
ylabel('Feeding level (with T=15C)')
subplot(3,1,3)
plot(R,f_l1,'b','LineWidth',2); hold on;
plot(R,f_l2,'r','LineWidth',2); hold on;
plot(R,f_l3,'k','LineWidth',2); hold on;
title('Large/Adults 630mm 2.5kg')
xlabel('prey density')
print('-dpng','feeding_level_T15.png')

%
T=25;
Acmax_s = exp(0.063*(T-15.0)) * h * M_s^(3/4);
Acmax_m = exp(0.063*(T-15.0)) * h * M_m^(3/4);
Acmax_l = exp(0.063*(T-15.0)) * h * M_l^(3/4);
beta_s1 = exp(0.063*T-15.0) * flev1 * M_s^(q);
beta_m1 = exp(0.063*T-15.0) * flev1 * M_m^(q);
beta_l1 = exp(0.063*T-15.0) * flev1 * M_l^(q);
beta_s2 = exp(0.063*T-15.0) * flev2 * M_s^(q);
beta_m2 = exp(0.063*T-15.0) * flev2 * M_m^(q);
beta_l2 = exp(0.063*T-15.0) * flev2 * M_l^(q);
beta_s3 = exp(0.063*T-15.0) * flev3 * M_s^(q);
beta_m3 = exp(0.063*T-15.0) * flev3 * M_m^(q);
beta_l3 = exp(0.063*T-15.0) * flev3 * M_l^(q);

f_s1 = beta_s1.*R ./ (Acmax_s + beta_s1.*R);
f_s2 = beta_s2.*R ./ (Acmax_s + beta_s2.*R);
f_s3 = beta_s3.*R ./ (Acmax_s + beta_s3.*R);
f_m1 = beta_m1.*R ./ (Acmax_m + beta_m1.*R);
f_m2 = beta_m2.*R ./ (Acmax_m + beta_m2.*R);
f_m3 = beta_m3.*R ./ (Acmax_m + beta_m3.*R);
f_l1 = beta_l1.*R ./ (Acmax_l + beta_l1.*R);
f_l2 = beta_l2.*R ./ (Acmax_l + beta_l2.*R);
f_l3 = beta_l3.*R ./ (Acmax_l + beta_l3.*R);

figure(6)
subplot(3,1,1)
plot(R,f_s1,'b','LineWidth',2); hold on;
plot(R,f_s2,'r','LineWidth',2); hold on;
plot(R,f_s3,'k','LineWidth',2); hold on;
title('T=25C')
legend('g=3500','g=4500','g=5500')
subplot(3,1,2)
plot(R,f_m1,'b','LineWidth',2); hold on;
plot(R,f_m2,'r','LineWidth',2); hold on;
plot(R,f_m3,'k','LineWidth',2); hold on;
ylabel('Feeding level')
subplot(3,1,3)
plot(R,f_l1,'b','LineWidth',2); hold on;
plot(R,f_l2,'r','LineWidth',2); hold on;
plot(R,f_l3,'k','LineWidth',2); hold on;
xlabel('prey density')









