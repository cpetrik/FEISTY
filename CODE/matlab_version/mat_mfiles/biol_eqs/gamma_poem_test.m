

clear all

K=1;
Z=500.0;
%nu=0.0012043900125042264;
d=4.854471554792261e-15;
B=2.619938604590845;
mu = d/B;
k=0.0;
if k==0.0
    kap=K;
else
    kap=1;
end
nu = 0:1;
gg1 = ((kap.*nu) - mu)./(1-(Z.^(1-(mu./(kap.*nu)))));

gg2 = ((kap.*nu) - mu)./(1-(Z.^((1-mu)./(kap.*nu))));

z = 1/Z;
gg3 = ((kap.*nu) - mu)./(1-(z.^(1-(mu./(kap.*nu)))));

%%
figure
%subplot(2,2,1)
plot(nu,gg1,'LineWidth',2); hold on;
ylabel('gamma = somatic growth')
xlabel('nu = energy for growth')
%print('-dpng','gamma_eq.png')

figure
plot(nu,gg1,'LineWidth',2); hold on;
plot(nu,gg3,'LineWidth',2); hold on;
legend('Z=1000','Z=1/1000')
ylabel('gamma = somatic growth')
xlabel('nu = energy for growth')
%print('-dpng','gamma_eq_fix.png')

%% DeRoos 2008

nu = 2;
nat = 0:1e-3:0.01;
pred = 0:2e-2:0.3;
g = NaN*ones(length(nat),length(pred));
for i=1:length(nat)
    for j=1:length(pred)
        mu = nat(i)+pred(j);
        Mu(i,j) = nat(i)+pred(j);
        gg3 = ((kap.*nu) - mu)./(1-(z.^(1-(mu./(kap.*nu)))));
        g(i,j) = (nu-mu) ./ (1 - z.^(1-mu./nu));
    end
end

[natG,predG] = meshgrid(nat,pred);
%%
figure
pcolor(natG,predG,g'); hold on;
shading flat
xlabel('natural mort')
ylabel('pred mort')
colorbar

figure
plot(Mu(:),g(:),'.'); hold on;
xlabel('total mort')
ylabel('gamma')





