% Look at output from orig Matlab code 
% to compare with new code

load('/Volumes/GFDL/CSV/advect_tests/POEM_adv_diff_test_Across_Arc_dt1hr_velH_monthly_b100.mat')

pdiff
m2=sum(TF2(:).*area(:));
m0=sum(TF0(:).*area(:));
100*(m2 - m0)/m0