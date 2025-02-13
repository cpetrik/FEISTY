% Megrey consumption vs. Ken Andersen's

k1=0.15;
k2=0.35;
k3=0.78;
k4=1.15;
h=60.0 / 365.0;
flev=0.5;
q=1;

L_s = 10^((log10(2)+log10(20))/2); 
L_m = 10^((log10(20)+log10(200))/2); 
L_l = 10^((log10(200)+log10(2000))/2); 

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

T = -2:0.2:37;

%% Without temp-dep part
R=0:0.1:10;
%Megrey et al.
Mcmax_s = 0.642 * M_s^(-0.256);
Mcmax_m = 0.642 * M_m^(-0.256);
Mcmax_l = 0.642 * M_l^(-0.256);

Mcon_s1 = (Mcmax_s .* 1.0 .* R ./ k1) ./ (1+ (1.0 .* R ./ k1));
Mcon_s2 = (Mcmax_s .* 0.1 .* R ./ k2) ./ (1+ (1.0 .* R ./ k1));
Mcon_m1 = (Mcmax_m .* 1.0 .* R ./ k1) ./ (1+ (1.0 .* R ./ k4));
Mcon_m2 = (Mcmax_m .* 0.5 .* R ./ k2) ./ (1+ (1.0 .* R ./ k4));
Mcon_l1 = (Mcmax_l .* 1.0 .* R ./ k1) ./ (1+ (1.0 .* R ./ k3));
Mcon_l2 = (Mcmax_l .* 0.5 .* R ./ k2) ./ (1+ (1.0 .* R ./ k3));

%Andersen
Acmax_s = h * M_s^(3/4);
Acmax_m = h * M_m^(3/4);
Acmax_l = h * M_l^(3/4);

beta_s = flev * M_s^(q);
beta_m = flev * M_m^(q);
beta_l = flev * M_l^(q);

Acon_s = Acmax_s .* (beta_s .* R) ./ (Acmax_s + beta_s.*R);
Acon_m = Acmax_m .* (beta_m .* R) ./ (Acmax_m + beta_s.*R);
Acon_l = Acmax_l .* (beta_l .* R) ./ (Acmax_l + beta_s.*R);

Acon_s2 = Acmax_s .* (1 .* R) ./ (Acmax_s + 1.*R);
Acon_m2 = Acmax_m .* (1 .* R) ./ (Acmax_m + 1.*R);
Acon_l2 = Acmax_l .* (1 .* R) ./ (Acmax_l + 1.*R);

%%
figure(1)
plot(R,Mcon_s1,'b','LineWidth',2); hold on;
plot(R,Mcon_s2,'--b','LineWidth',2); hold on;
plot(R,Acon_s,'c','LineWidth',2); hold on;
plot(R,Acon_s2,'--c','LineWidth',2); hold on;
ylim([0 3])

figure(2)
plot(R,Mcon_m1,'r','LineWidth',2); hold on;
plot(R,Mcon_m2,'--r','LineWidth',2); hold on;
plot(R,Acon_m,'m','LineWidth',2); hold on;
plot(R,Acon_m2,'--m','LineWidth',2); hold on;
ylim([0 5])

figure(3)
plot(R,Mcon_l1,'k','LineWidth',2); hold on;
plot(R,Mcon_l2,'--k','LineWidth',2); hold on;
plot(R,Acon_l,'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(R,Acon_l2,'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
ylim([0 5])




