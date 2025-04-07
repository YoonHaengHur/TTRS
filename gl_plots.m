%% gl_fixed_d
clear
load('gl_fixed_d.mat')
log_errors = log10(errors)';
log_errors_mean = mean(log_errors);
log_errors_std = std(log_errors);

center = [mean(log10(NN)), mean(log_errors_mean)];
errorbar(log10(NN), log_errors_mean, log_errors_std, '.-b')
hold on
plot(log10(NN), -0.5*(log10(NN)-center(1))+center(2), '--r')
xlabel('log10(N)','FontSize',14)
ylabel('log10(relative error)','FontSize',14)
xlim([log10(NN(1)), inf])
legend('log10(error)','log10($1/\sqrt{N}$)','Interpreter','latex')
hold off
set(gca,'box','off')
h=gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'gl_fixed_d','-dpdf','-r0')
close

%% gl_fixed_N
clear
load('gl_fixed_N.mat')
errors = errors';
errors_mean = mean(errors);
errors_std = std(errors);

errorbar(dd, errors_mean, errors_std, '.-b')
xlabel('d','FontSize',14)
ylabel('relative error','FontSize',14)
xlim([dd(1), inf])
hold off
set(gca,'box','off')
h=gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'gl_fixed_N','-dpdf','-r0')
close
