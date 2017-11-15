A = [80	  7.9368
160	 29.8970
240	 65.4789
320	116.0674
400	177.6208
480	253.6504
560	344.2471
642	448.6176];

x = A(:,1);
y = A(:,2);
slope = 1.5*10^(-3)*x.^2;

figure
loglog(x,slope,'-m','LineWidth',1);
hold on;
loglog(x,y,'-ob','LineWidth',2);
hold off;
xlim([60 800]);
ylim([6 800]);
grid on;
set(gca,'FontSize',14)
xlabel('number of atoms','fontsize',18)
ylabel('execution time (sec.)','fontsize',18)
legend({'slope 2'},'fontsize',14,'location','northwest')
saveas(gcf,'plot.png')

