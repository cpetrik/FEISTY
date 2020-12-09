function plotPoemf(param, result)

y = result.y;
R = result.R;
B = result.B;
t = result.t;

xlimit = [min(param.wc(param.ixFish))/10 max(param.wu)];
Bin = floor(0.8*length(y));
yend = mean(y(Bin:end,:));
[f, mortpred] = calcEncounter(yend', param);
wc = param.wc;

clf
figure(1)
subplot(3,1,1)
semilogx(param.wc(param.ix1(1):param.ix2(1)),yend(param.ix1(1):param.ix2(1)),'linewidth',2)
hold on
for ii = 2 : param.nSpecies
semilogx(param.wc(param.ix1(ii):param.ix2(ii)),yend(param.ix1(ii):param.ix2(ii)),'linewidth',2)
end
semilogx(param.wc(param.ix1(1):param.ix2(1)),yend(param.ix1(1):param.ix2(1)),'Color',[0  0.4470  0.7410],'linewidth',2)
xlabel('central weight (grams)')
ylabel('Biomass')
xlim(xlimit)
ylim([1e-6 5]);
legend('Spel','Meso','Lpel','Bpel', 'Dem')
hold off

subplot(3,1,2)
semilogx(param.wc(param.ix1(1):param.ix2(1)),f(end,param.ix1(1):param.ix2(1)),'linewidth',2)
hold on
for ii = 2 : param.nSpecies
semilogx(param.wc(param.ix1(ii):param.ix2(ii)),f(end,param.ix1(ii):param.ix2(ii)),'linewidth',2)
end
semilogx(param.wc(param.ix1(1):param.ix2(1)),f(end,param.ix1(1):param.ix2(1)),'Color',[0  0.4470  0.7410],'linewidth',2)
semilogx(xlimit, 0.2*[1 1], 'k--')
ylim([0 1])
xlim(xlimit)
xlabel('central weight (grams)')
ylabel('Feeding lvl.')
hold off

subplot(3,1,3)
semilogx(param.wc(param.ix1(1):param.ix2(1)),mortpred(end,param.ix1(1):param.ix2(1)),'linewidth',2)
hold on
for ii = 2 : param.nSpecies
semilogx(param.wc(param.ix1(ii):param.ix2(ii)),mortpred(end,param.ix1(ii):param.ix2(ii)),'linewidth',2)
end
semilogx(param.wc(param.ix1(1):param.ix2(1)),mortpred(end,param.ix1(1):param.ix2(1)),'Color',[0  0.4470  0.7410],'linewidth',2)
ylim([0 1.2*max(mortpred)])
xlim(xlimit)
ylabel('Predation mort.')
xlabel('central weight (grams)')
hold off
