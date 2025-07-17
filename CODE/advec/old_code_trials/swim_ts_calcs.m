L_s = 10.0; % small
L_m = 200.0; % medium
L_l = 1.0e3;% large

%%! Mass from length using Andersen & Beyer 2013
% Convert from mm to cm and use their const coeff = 0.01g/cm3
M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

wgt = M_m; 	
T = 30.0;
w = ((3.9*wgt.^0.13 * exp(0.149*T)) /100);
dw = w * 60 * 60 *24

wgt2 = M_l; 	
w2 = ((3.9*wgt2.^0.13 * exp(0.149*T)) /100);
dw2 = w2 * 60 * 60 *24

min(dyte(:))
mean(dyte(:))
median(dyte(:))

min(dxtn(:))
mean(dxtn(:))
median(dxtn(:))

min(sqrt(dat(:)))
mean(sqrt(dat(:)))
median(sqrt(dat(:)))

dist = min(sqrt(dat(:)));
dist / (w * 60 * 60)
dist / (w2 * 60 * 60)