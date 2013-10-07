t1 = [0:150];
pot = exp(-t1/17);

t2 = [-150:0];
dep = -1*exp(t2/34);

figure
plot([-150 150],[0 0],'LineStyle','--','color',[0.7 0.7 0.7])
hold on
plot([0 0],[-1 1],'LineStyle','--','color',[0.7 0.7 0.7])
plot([0:150],pot,'LineWidth',2,'color','k')
plot([-150:0],dep,'LineWidth',2,'color','k')
ylabel('Synaptic changes')
xlabel('Spike timing [msec]')