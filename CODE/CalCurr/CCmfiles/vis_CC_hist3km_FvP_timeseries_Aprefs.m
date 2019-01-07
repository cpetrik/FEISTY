% Visualize output of POEM
% Historic 1988-2010 with 3km model
% Time series plots of F vs. P with different A prefs

clear all
close all

% Fish data
ppath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/Pref_tests/CC/';
%harv = 'pristine';
%tharv = 'F=0';
harv = 'All_fish01';
tharv = 'All fish F=0.1';
cpath = '/Volumes/GFDL/NC/Matlab_new_size/bio_rates/CCE/';

%% Colors
cm4=[0 0 1;...    %b
    1 0 0;...     %r
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black
    
set(groot,'defaultAxesColorOrder',cm4);

%%
ad=0.1:0.1:1;
nt=12*23;
all_biom = NaN*ones(length(ad),nt,4);
for n=1:length(ad)
    A = ad(n);
    ta = num2str(1000+int64(100*A));
    cfile = ['Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A',ta(2:end),...
        '_Sm025_nmort1_BE08_noCC_RE00100'];
    fpath = ['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/CalCurr/'];
    load([fpath 'Means_Historic3km_' harv '_' cfile '.mat']);
    
    
    %% Plots in time
    t = 1:length(sp_tmean); %time;
    y = 1988 + (t-1)/12;
    
    F = sf_tmean+mf_tmean;
    P = sp_tmean+mp_tmean+lp_tmean;
    D = sd_tmean+md_tmean+ld_tmean;
    B = b_tmean;
    
    all_biom(n,:,1) = F;
    all_biom(n,:,2) = P;
    all_biom(n,:,3) = D;
    all_biom(n,:,4) = B;
    
    %%
    figure(1)
    clf
    yyaxis left
    plot(y,(P),'Linewidth',2);
    ylabel('P Biomass (g m^-^2)')
    xlim([y(1) y(end)])
    ylim([0 25])
    yyaxis right
    plot(y,(F),'Linewidth',2);
    ylabel('F Biomass (g m^-^2)')
    xlim([y(1) y(end)])
    ylim([0 25])
    xlabel('Time (y)')
    title(['Historic 3km ' tharv])
    stamp(['A = ' num2str(A)])
    print('-dpng',[ppath 'Historic3km_',harv,'_',cfile,'_FvP.png'])
    
    if (n<10)
        f2=figure(2)
        subplot(3,3,n)
        yyaxis left
        plot(y,(P),'Linewidth',2);
        ylabel('P Biomass (g m^-^2)')
        xlim([y(1) y(end)])
        ylim([0 25])
        yyaxis right
        plot(y,(F),'Linewidth',2);
        ylabel('F Biomass (g m^-^2)')
        xlim([y(1) y(end)])
        ylim([0 25])
        xlabel('Time (y)')
        title(['A = ' num2str(A)])
        if (n==9)
            stamp(harv)
        end
    end
    
end
print(f2,'-dpng',[ppath 'Historic3km_',harv,'_Aprefs_FvP.png'])

save([cpath 'Historic3km_',harv,'_Aprefs_AllBiom.png'],'all_biom','ad','cfile')
