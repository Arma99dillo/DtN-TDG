function PlotConvergence(ni,nf,L2Error,H1Error,ticks)

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


figure()
semilogy(ni:nf,L2Error,'.-r',ni:nf,H1Error,'*--b'); grid
xlim([ni nf]); xticks(ticks);
LL = legend('$L^2$ error','$H^1$ error','FontSize', 14);
set(LL, 'Interpreter', 'latex');
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'TickLabelInterpreter', 'latex')
a = get(gca,'YTickLabel');  
set(gca,'YTickLabel',a,'fontsize',14,'TickLabelInterpreter', 'latex')
xlabel('Number of local plane wave basis functions','FontSize',18, 'Interpreter','latex')
ylabel('Error','FontSize',18, 'Interpreter','latex')

end