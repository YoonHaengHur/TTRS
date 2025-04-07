%% ising_fixed_d
clear
load('ising_fixed_d.mat')
log_errors1 = log10(errors1)';
log_errors2 = log10(errors2)';

h = figure;
shadedErrorBar(log10(NN), log_errors1, {@mean,@std}, 'lineprops', '+-r','patchSaturation',0.075)
shadedErrorBar(log10(NN), log_errors2, {@mean,@std}, 'lineprops', 'x--b','patchSaturation',0.075)
xlim([log10(NN(1)), inf])
legend('sketch1','sketch2')
xlabel('log10(N)','FontSize',14)
ylabel('log10(relative error)','FontSize',14)

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'ising_fixed_d','-dpdf','-r0')
close

%% ising5_fixed_d
clear
load('ising5_fixed_d.mat')
log_errors1 = log10(errors1)';
log_errors2 = log10(errors2)';

h = figure;
shadedErrorBar(log10(NN), log_errors1, {@mean,@std}, 'lineprops', '+-r','patchSaturation',0.075)
shadedErrorBar(log10(NN), log_errors2, {@mean,@std}, 'lineprops', 'x--b','patchSaturation',0.075)
xlim([log10(NN(1)), inf])
legend('sketch1','sketch2')
xlabel('log10(N)','FontSize',14)
ylabel('log10(relative error)','FontSize',14)
hold off

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'ising5_fixed_d','-dpdf','-r0')
close

%% ising_fixed_N
clear
load('ising_fixed_N.mat')
errors1 = errors1';
errors2 = errors2';

h = figure;
shadedErrorBar(dd, errors1, {@mean,@std}, 'lineprops', '+-r','patchSaturation',0.075)
shadedErrorBar(dd, errors2, {@mean,@std}, 'lineprops', 'x--b','patchSaturation',0.075)
xlim([dd(1), inf])
legend('sketch1','sketch2')
xlabel('d','FontSize',14)
ylabel('relative error','FontSize',14)
hold off

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'ising_fixed_N','-dpdf','-r0')
close

%% ising5_fixed_N
clear
load('ising5_fixed_N.mat')
errors1 = errors1';
errors2 = errors2';

h = figure;
shadedErrorBar(dd, errors1, {@mean,@std}, 'lineprops', '+-r','patchSaturation',0.075)
shadedErrorBar(dd, errors2, {@mean,@std}, 'lineprops', 'x--b','patchSaturation',0.075)
xlim([dd(1), inf])
legend('sketch1','sketch2')
xlabel('d','FontSize',14)
ylabel('relative error','FontSize',14)
hold off

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'ising5_fixed_N','-dpdf','-r0')
close

