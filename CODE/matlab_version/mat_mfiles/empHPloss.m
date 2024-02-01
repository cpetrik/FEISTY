% HP loss is an empirical fitted fn of biomass and temp

% USE HISTORIC BECAUSE COVERS A GREATER RANGE OF TEMPERATURES

%% Medium zoo
% Historic
%linear, log
% lMZh <- lm(Lmzl ~ Lmzb + pt, data=HZ)
% Coefficients:
%               Estimate Std. Error t value Pr(>|t|)    
% (Intercept) -2.617e+00  3.291e-04   -7951   <2e-16 ***
% Lmzb         1.989e+00  5.162e-04    3854   <2e-16 ***
% pt           1.732e-02  9.567e-06    1810   <2e-16 ***
% ---
% Adjusted R-squared:  0.9446
dZmd = 10 .^ (-2.617 + 1.989.*log10(Zmd+eps) + 1.732e-2.*Tp);

% Climatol
%linear, log
% lMZh <- lm(Lmzl ~ Lmzb + pt, data=HZ)
% Coefficients:
%               Estimate Std. Error t value Pr(>|t|)    
% (Intercept) -2.616e+00  1.957e-03 -1337.2   <2e-16 ***
% Lmzb         2.035e+00  3.036e-03   670.4   <2e-16 ***
% pt           1.646e-02  7.429e-05   221.5   <2e-16 ***
% ---
% Adjusted R-squared:  0.932
dZmd = 10 .^ (-2.616 + 2.035.*log10(Zmd+eps) + 1.646e-2.*Tp);

%% Large zoo
% Historic
%linear, log
% lLZh <- lm(Llzl ~ Llzb + pt, data=HZ)
% summary(lLZh) 
% Coefficients:
%               Estimate Std. Error t value Pr(>|t|)    
% (Intercept) -2.954e+00  2.627e-04  -11244   <2e-16 ***
% Llzb         2.228e+00  4.294e-04    5190   <2e-16 ***
% pt           2.556e-02  1.285e-05    1990   <2e-16 ***
% ---
% Adjusted R-squared:  0.9673 
dZlg = 10 .^ (-2.954 + 2.228.*log10(Zlg+eps) + 2.556e-2.*Tp);

% Climatol
%linear, log
% lLZh <- lm(Llzl ~ Llzb + pt, data=HZ)
% summary(lLZh) 
% Coefficients:
%               Estimate Std. Error t value Pr(>|t|)    
% (Intercept) -2.940e+00  1.279e-03 -2298.5   <2e-16 ***
% Llzb         2.210e+00  1.977e-03  1118.1   <2e-16 ***
% pt           2.599e-02  6.437e-05   403.8   <2e-16 ***
% ---
% Adjusted R-squared:  0.9737
dZlg = 10 .^ (-2.940 + 2.210.*log10(Zlg+eps) + 2.599e-2.*Tp);

%% One mesozoo
% Historic
%linear, log
dZmeso = 10 .^ (-2.925 + 1.964.*log10(Zmeso+eps) + 1.958e-2.*Tp);
    

