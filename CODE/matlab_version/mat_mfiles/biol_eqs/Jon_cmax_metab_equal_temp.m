% Test for Jonathan's O2 model
% What temp does orig total metab = 60% Cmax?

clear 
close all

prey=0:0.005:0.3;
temp=15;

%!Individual Mass (g)
M_s = 10^((log10(0.001)+log10(0.5))/2);
M_m = 10^((log10(0.5)+log10(250))/2);
M_l = 10^((log10(250)+log10(125000))/2);

%! Body lengths (mm)
% Convert from mm to cm and use their const coeff = 0.01g/cm3
L_s = 10.0 * (M_s/0.01)^(1/3); % small
L_m = 10.0 * (M_m/0.01)^(1/3); % medium
L_l = 10.0 * (M_l/0.01)^(1/3); % large

m=[M_s; M_m; M_l];

%Me
Ccmax = (exp(0.063*(temp-10.0)) .* 20 .* m.^(-0.25)) ./365.0;

Cmet = 0.2 * (exp(0.08555*(temp-10.0)) .* 20 .* m.^(-0.175)) ./365.0;

% TEMP-DEP
mass=[M_s; M_m; M_l];
temp = -2:40;
temp2 = temp+273;

for s=1:3
    m = mass(s);
    
   
    
    Ccmax = (exp(0.063*(temp-10.0)) .* 20 .* m.^(-0.25)) ./365.0;
    
    
    Cmet = 0.2 * (exp(0.0905*(temp-10.0)) .* 20 .* m.^(-0.175)) ./365.0;
    
    
    f5=figure(5);
     
    subplot(3,1,s)
    plot(temp,log10(Cmet),'-b','MarkerSize',15,'LineWidth',2); hold on;
     plot(temp,log10(0.7*Ccmax),'-','color',[0 0.5 0.75],'MarkerSize',15,'LineWidth',2); hold on;
    xlim([-2 35])
    if (s==1)
%         legend('K&H','Hart','J&C','mizer','me')
%         legend('location','northwest')
        
    end
    
    
    
end
