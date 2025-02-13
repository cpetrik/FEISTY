A = 4.39;
fc = 0.2;
f0 = 0.6;
epsassim = 0.6;
n = 3/4;
 
w = logspace(-3, 5);
 
AvailEnergy = A*w.^n;
Consumption = A / (epsassim*(f0-fc)) * w.^n;
 
 
loglog(w, AvailEnergy, w, Consumption)
xlabel('Weight (g)')
ylabel('Rate (g/yr)')
legend({'Available energy','Consumption'})