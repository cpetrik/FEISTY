function plotPoemf(param, result)

y = result.y;
R = result.R;
B = result.B;
t = result.t;

xlimit = [min(param.wc(param.ixFish))/10 max(param.wu)];
Bin = floor(0.8*length(y));
yend = mean(y(Bin:end,:));

% Use model result instead of calc
%[f, mortpred] = calcEncounter(yend', param);
f = result.f;
mortpred = result.mortpred;
wc = param.wc;

% colors
cm6=[1 0 0;...    %r
    0.5 0 0;...   %maroon
    0 0 0.75;...  %b
    0 0.5 0.75;... %med blue
    0 0.7 0;...   %g
    0 0 0];       %black

set(groot,'defaultAxesColorOrder',cm6);

%% Why are mid-water fish feeding level and mortality > 0?
%Change so if not in the model (or biomass = 0), then 0
Bend = yend;
Bend(Bend > 0) = 1;
f = result.f .* Bend;
mortpred = result.mortpred .* Bend;

%%
clf
figure(1)
subplot(3,1,1)
%semilogx(param.wc(param.ix1(1):param.ix2(1)),yend(param.ix1(1):param.ix2(1)),'linewidth',2)
%hold on
for ii = 1 : param.nSpecies
    semilogx(param.wc(param.ix1(ii):param.ix2(ii)),yend(param.ix1(ii):param.ix2(ii)),'linewidth',2); hold on
end
%semilogx(param.wc(param.ix1(1):param.ix2(1)),yend(param.ix1(1):param.ix2(1)),'Color',[0  0.4470  0.7410],'linewidth',2)
xlabel('central weight (grams)')
ylabel('Biomass')
xlim(xlimit)
ylim([1e-6 6]);
hold off

subplot(3,1,2)
% semilogx(param.wc(param.ix1(1):param.ix2(1)),f(end,param.ix1(1):param.ix2(1)),'linewidth',2)
% hold on
for ii = 1 : param.nSpecies
    semilogx(param.wc(param.ix1(ii):param.ix2(ii)),f(end,param.ix1(ii):param.ix2(ii)),'linewidth',2); hold on;
end
%semilogx(param.wc(param.ix1(1):param.ix2(1)),f(end,param.ix1(1):param.ix2(1)),'Color',[0  0.4470  0.7410],'linewidth',2)
semilogx(xlimit, 0.2*[1 1], 'k--')
ylim([0 1])
xlim(xlimit)
xlabel('central weight (grams)')
ylabel('Feeding lvl.')
hold off

subplot(3,1,3)
% semilogx(param.wc(param.ix1(1):param.ix2(1)),mortpred(end,param.ix1(1):param.ix2(1)),'linewidth',2)
% hold on
for ii = 1 : param.nSpecies
    semilogx(param.wc(param.ix1(ii):param.ix2(ii)),mortpred(end,param.ix1(ii):param.ix2(ii)),'linewidth',2); hold on;
end
% semilogx(param.wc(param.ix1(1):param.ix2(1)),mortpred(end,param.ix1(1):param.ix2(1)),'Color',[0  0.4470  0.7410],'linewidth',2)
ylim([0 1.2*max(mortpred)])
xlim(xlimit)
ylabel('Predation mort.')
xlabel('central weight (grams)')
legend('SmPel','MesPel','LgPel','BathPel', 'Dem')
hold off

%%
ltext = {'SmPel','MesPel','LgPel','BathPel', 'Dem'};
figure(2)
clf
for ii = 1 : param.nSpecies
    subplot(3,2,ii)
    semilogx(param.wc(param.ix1(ii):param.ix2(ii)),yend(param.ix1(ii):param.ix2(ii)),...
        'linewidth',2,'color',cm6(ii,:));
    xlabel('central weight (g)')
    ylabel('Biomass')
    xlim(xlimit)
    %ylim([1e-6 6]);
    title(ltext(ii))
end

figure(3)
clf
for ii = 1 : param.nSpecies
    subplot(3,2,ii)
    semilogx(param.wc(param.ix1(ii):param.ix2(ii)),f(end,param.ix1(ii):param.ix2(ii)),...
        'linewidth',2,'color',cm6(ii,:)); hold on
    semilogx(xlimit, 0.2*[1 1], 'k--')
    ylim([0 1])
    xlim(xlimit)
    xlabel('central weight (g)')
    ylabel('Feeding level')
    title(ltext(ii))
end
hold off

figure(4)
clf
for ii = 1 : param.nSpecies
    subplot(3,2,ii)
    semilogx(param.wc(param.ix1(ii):param.ix2(ii)),mortpred(end,param.ix1(ii):param.ix2(ii)),...
        'linewidth',2,'color',cm6(ii,:)); hold on;
    ylim([0 1.1*max(mortpred)])
    xlim(xlimit)
    ylabel('Predation mort')
    xlabel('central weight (g)')
    title(ltext(ii))
end
hold off

