% Old way

M_s = 10^((log10(0.001)+log10(0.5))/2)          %0.0224
M_m = 10^((log10(0.5)+log10(250))/2)            %11.1803
M_l = 10^((log10(250)+log10(125000))/2)         %5590.2

%! Body lengths (mm)
% Convert from mm to cm and use their const coeff = 0.01g/cm3
L_s = 10.0 * (M_s/0.01)^(1/3) % 13.0766
L_m = 10.0 * (M_m/0.01)^(1/3) % 103.7891
L_l = 10.0 * (M_l/0.01)^(1/3) % 823.7745


%Correct way
M_s = (0.001-0.5)/(log10(0.001)-log10(0.5))     %0.1849 
M_m = (0.5-250)/(log10(0.5)-log10(250))         %92.4427
M_l = (250-125000)/(log10(250)-log10(125000))   %46221

%! Body lengths (mm)
% Convert from mm to cm and use their const coeff = 0.01g/cm3
L_s = 10.0 * (M_s/0.01)^(1/3) % 26.4424
L_m = 10.0 * (M_m/0.01)^(1/3) % 209.8734
L_l = 10.0 * (M_l/0.01)^(1/3) % 1.6658e+03