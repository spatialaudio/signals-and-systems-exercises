clear all;
close all;
clc;
set(0,'defaulttextinterpreter','latex');
figure1 = figure('PaperUnits','centimeters','Units','normalized');

%matplotlib default colors
C0 = sscanf('1f77b4','%2x%2x%2x',[1 3])/255;
C1 = sscanf('ff7f0e','%2x%2x%2x',[1 3])/255;
C2 = sscanf('2ca02c','%2x%2x%2x',[1 3])/255;
C3 = sscanf('d62728','%2x%2x%2x',[1 3])/255;
C7 = sscanf('7f7f7f','%2x%2x%2x',[1 3])/255;

LW = 3;
N = 10000;
n=0:N-1;

subplot(2,1,1)
x1 = sin(2*pi*2/N*n);
plot(n,x1,'LineWidth',LW, 'color', C0), hold on
plot(n+1000,x1,'LineWidth',LW, 'color', C1)
legend('x(t)','y(t)')
axis([5000 N-1 -1 1])
set(gca,'XTick',0)
xlabel('$t \rightarrow $')
annotation(figure1,'doublearrow',[0.32 0.48],[0.9 0.9], 'LineWidth',4, 'color',C3);
title('Phasenlaufzeit (zeitlicher Versatz der Signale)')
text(6700,0.5,'$\tau_\mathrm{PD}$','FontSize',18, 'color', C3)

subplot(2,1,2)
a=1/200;
Mod = a/sqrt(pi)*exp(-(a*[-N/2:1:N/2-1]).^2);
Mod = Mod/max(abs(Mod));
x2 = sin(2*pi*100/N*n).*Mod;
plot(n,+Mod,'LineWidth',LW, 'color',C0), hold on
plot(n,-Mod,'LineWidth',LW, 'color', C0),
plot(n,x2,'LineWidth',LW, 'color', C0)
plot(n+500,+Mod,'LineWidth',LW, 'color', C1),
plot(n+500,-Mod,'LineWidth',LW, 'color', C1),
plot(n+500,x2,'LineWidth',LW, 'color', C1),
hold off
axis([4000 6500 -1 1])
set(gca,'XTick',0)
xlabel('$t \rightarrow $')
annotation(figure1,'doublearrow',[0.435 0.595],[0.42 0.42], 'LineWidth',4, 'color',C2);
text(5200,0.7,'$\tau_\mathrm{GD}$','FontSize',18, 'color',C2)
title(['Gruppenlaufzeit (zeitlicher Versatz der Einh\"ullenden)'])

print -depsc2 PhasenGruppenlaufzeit.eps
