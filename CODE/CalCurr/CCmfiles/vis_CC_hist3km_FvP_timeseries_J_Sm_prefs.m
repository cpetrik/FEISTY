% Visualize output of POEM
% Historic 1988-2010 with 3km model
% Time series plots and pcolor of F vs. P with different J and Sm prefs

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
juv=0.1:0.1:1;
sml=0.1:0.1:1;
nt=12*23;

allF = NaN*ones(nt,length(sml),length(juv));
allP = NaN*ones(nt,length(sml),length(juv));
allD = NaN*ones(nt,length(sml),length(juv));
allB = NaN*ones(nt,length(sml),length(juv));

for n=1:7%:length(juv)
    J = juv(n);
    tj = num2str(1000+int64(100*J));
    for m=1:length(sml)
        Sm = sml(m);
        ts = num2str(1000+int64(100*Sm));
        
        cfile = ['Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J',tj(2:end),...
            '_A100_Sm',ts(2:end),'_nmort1_BE08_noCC_RE00100'];
        fpath = ['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/CalCurr/'];
        load([fpath 'Historic3km_' harv '_' cfile '_Means.mat']);
        
        
        %% Plots in time
        t = 1:length(sp_tmean); %time;
        y = 1988 + (t-1)/12;
        
        F = sf_tmean+mf_tmean;
        P = sp_tmean+mp_tmean+lp_tmean;
        D = sd_tmean+md_tmean+ld_tmean;
        B = b_tmean;
        
        allF(:,m,n) = F';
        allP(:,m,n) = P';
        allD(:,m,n) = D';
        allB(:,m,n) = B';
        
        if (n<4)
            if (m<10)
                f1=figure(1);
                subplot(3,3,m)
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
                title(['Sm = ' num2str(Sm)])
                if (m==9)
                    stamp(['J = ' num2str(J) ' ' harv])
                end
            end
        end
        
    end
    if (n<4)
        print(f1,'-dpng',[ppath 'Historic3km_',harv,'_J',tj(2:end),...
            '_SmPrefs_FvP.png'])
    end
end


%%
[sgrid,jgrid]=meshgrid([sml sml(end)+0.1],[juv juv(end)+0.1]);
nj=length(juv);
ns=length(sml);
mF = squeeze(nanmean(allF));
mF(ns+1,:) = nan;
mF(:,nj+1) = nan;
mP = squeeze(nanmean(allP));
mP(ns+1,:) = nan;
mP(:,nj+1) = nan;

fracPF = (mP) ./ (mP+mF);

%% Sum mean biom over stages
f2=figure(2);
subplot(2,2,1)
pcolor(jgrid,sgrid,log10(mF))
colorbar
caxis([-1 0])
title('log10 Mean F Biom (g m^-^2)')
xlabel('Sm pref')
ylabel('J pref')

subplot(2,2,2)
pcolor(jgrid,sgrid,log10(mP))
colorbar
caxis([0.5 1.5])
title('log10 Mean P Biom (g m^-^2)')
xlabel('Sm pref')
ylabel('J pref')

subplot(2,2,3)
pcolor(jgrid,sgrid,fracPF)
colorbar
caxis([0.9 1])
title('Fraction Pelagics')
xlabel('Sm pref')
ylabel('J pref')
stamp(harv)
print(f2,'-dpng',[ppath 'Historic3km_',harv,'_J_Sm_prefs_FvP.png'])


%%
save([cpath 'Historic3km_',harv,'_J_Sm_prefs_AllBiom.png'],'juv',...
    'sml','cfile','allF','allP','allD','allB')


