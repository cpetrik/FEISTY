%% Find negative values
bio=NaN*ones(9,365);
con=NaN*ones(8,365);
die=NaN*ones(8,365);
gamma=NaN*ones(8,365);
rec=NaN*ones(8,365);

for n=1:365
    %bio
    id = (S_Sml_p(:,n) < 0);
    bio(1,n) = sum(id);
    id = (S_Sml_f(:,n) < 0);
    bio(2,n) = sum(id);
    id = (S_Sml_d(:,n) < 0);
    bio(3,n) = sum(id);
    id = (S_Med_p(:,n) < 0);
    bio(4,n) = sum(id);
    id = (S_Med_f(:,n) < 0);
    bio(5,n) = sum(id);
    id = (S_Med_d(:,n) < 0);
    bio(6,n) = sum(id);
    id = (S_Lrg_p(:,n) < 0);
    bio(7,n) = sum(id);
    id = (S_Lrg_d(:,n) < 0);
    bio(8,n) = sum(id);
    id = (S_Bent_bio(:,n) < 0);
    bio(9,n) = sum(id);
    %con
    id = (S_Sml_p_con(:,n) < 0);
    con(1,n) = sum(id);
    id = (S_Sml_f_con(:,n) < 0);
    con(2,n) = sum(id);
    id = (S_Sml_d_con(:,n) < 0);
    con(3,n) = sum(id);
    id = (S_Med_p_con(:,n) < 0);
    con(4,n) = sum(id);
    id = (S_Med_f_con(:,n) < 0);
    con(5,n) = sum(id);
    id = (S_Med_d_con(:,n) < 0);
    con(6,n) = sum(id);
    id = (S_Lrg_p_con(:,n) < 0);
    con(7,n) = sum(id);
    id = (S_Lrg_d_con(:,n) < 0);
    con(8,n) = sum(id);
    %die
    id = (S_Sml_p_die(:,n) < 0);
    die(1,n) = sum(id);
    id = (S_Sml_f_die(:,n) < 0);
    die(2,n) = sum(id);
    id = (S_Sml_d_die(:,n) < 0);
    die(3,n) = sum(id);
    id = (S_Med_p_die(:,n) < 0);
    die(4,n) = sum(id);
    id = (S_Med_f_die(:,n) < 0);
    die(5,n) = sum(id);
    id = (S_Med_d_die(:,n) < 0);
    die(6,n) = sum(id);
    id = (S_Lrg_p_die(:,n) < 0);
    die(7,n) = sum(id);
    id = (S_Lrg_d_die(:,n) < 0);
    die(8,n) = sum(id);
    %gamma/rep
    id = (S_Sml_p_gamma(:,n) < 0);
    gamma(1,n) = sum(id);
    id = (S_Sml_f_gamma(:,n) < 0);
    gamma(2,n) = sum(id);
    id = (S_Sml_d_gamma(:,n) < 0);
    gamma(3,n) = sum(id);
    id = (S_Med_p_gamma(:,n) < 0);
    gamma(4,n) = sum(id);
    id = (S_Med_f_rep(:,n) < 0);
    gamma(5,n) = sum(id);
    id = (S_Med_d_gamma(:,n) < 0);
    gamma(6,n) = sum(id);
    id = (S_Lrg_p_rep(:,n) < 0);
    gamma(7,n) = sum(id);
    id = (S_Lrg_d_rep(:,n) < 0);
    gamma(8,n) = sum(id);
    %rec
    id = (S_Sml_p_rec(:,n) < 0);
    rec(1,n) = sum(id);
    id = (S_Sml_f_rec(:,n) < 0);
    rec(2,n) = sum(id);
    id = (S_Sml_d_rec(:,n) < 0);
    rec(3,n) = sum(id);
    id = (S_Med_p_rec(:,n) < 0);
    rec(4,n) = sum(id);
    id = (S_Med_f_rec(:,n) < 0);
    rec(5,n) = sum(id);
    id = (S_Med_d_rec(:,n) < 0);
    rec(6,n) = sum(id);
    id = (S_Lrg_p_rec(:,n) < 0);
    rec(7,n) = sum(id);
    id = (S_Lrg_d_rec(:,n) < 0);
    rec(8,n) = sum(id);
    
end
%%
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Mat_runs/';
cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D100_nmort0_BE05_CC275_RE0500';
ppath = [pp cfile '/'];

figure
plot(1:365,bio)
ylim([0 10])
legend('SP','MP','LP','SF','MF','SD','MD','LD','B')
legend('location','northwest')
print('-dpng',[ppath 'biotest.png'])

figure
plot(1:365,con)
ylim([0 10])
legend('SP','MP','LP','SF','MF','SD','MD','LD')
legend('location','northwest')
print('-dpng',[ppath 'contest.png'])

figure
plot(1:365,die)
ylim([0 10])
legend('SP','MP','LP','SF','MF','SD','MD','LD')
legend('location','northwest')
print('-dpng',[ppath 'dietest.png'])

figure
plot(1:365,gamma)
ylim([0 10])
legend('SP','MP','LP','SF','MF','SD','MD','LD')
legend('location','northwest')
print('-dpng',[ppath 'gammatest.png'])

figure
plot(1:365,rec)
ylim([0 10])
legend('SP','MP','LP','SF','MF','SD','MD','LD')
legend('location','northwest')
print('-dpng',[ppath 'rectest.png'])
