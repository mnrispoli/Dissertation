clear all
close all
clc

prefix=''
Dat8s=csvread([prefix 'x_norm_8s_dat.csv'])
Thry8s=csvread([prefix 'x_norm_8s_thry.csv'])
Thry8sU0=csvread([prefix 'x_norm_U0_8s_thry.csv'])
Dat12s=csvread([prefix 'x_norm_12s_dat.csv'])
Thry12s=csvread([prefix 'x_norm_12s_thry.csv'])
Thry12sU0=csvread([prefix 'x_norm_U0_12s_thry.csv'])



figure(1)
subplot(1,2,1)
hold off
ps0=plot(Thry8sU0(:,1),Thry8sU0(:,2),'k--', ...
    'LineWidth',2)
hold on
ps8=plot(Thry8s(:,1),Thry8s(:,2),'-', ...
    'LineWidth',2)
pe8=errorbar(Dat8s(:,1),Dat8s(:,2),Dat8s(:,3),'o', ...
    'LineWidth',2, ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', 'White', ...
    'CapSize',0 ...
    )

pe8.Color=ps8.Color;
ps0.Color=ps8.Color;

apar=0.6;
ps8.Color=(apar).*(ps8.Color)+([1 1 1].*(1-apar))
ps0.Color=(apar).*(ps0.Color)+([1 1 1].*(1-apar))


set(gcf,'color','white')
grid on
xlabel('Disorder W(J)')
legend('U=0','U=2.7J','Data')
ylabel('\Delta x (sites)')
title('L=8 Comparison')


subplot(1,2,2)
hold off
ps0=plot(Thry12sU0(:,1),Thry12sU0(:,2),'k--', ...
    'LineWidth',2)
hold on
ps8=plot(Thry12s(:,1),Thry12s(:,2),'-', ...
    'LineWidth',2)
pe8=errorbar(Dat12s(:,1),Dat12s(:,2),Dat12s(:,3),'o', ...
    'LineWidth',2, ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', 'White', ...
    'CapSize',0 ...
    )

pe8.Color=ps8.Color;
ps0.Color=ps8.Color;

apar=0.6;
ps8.Color=(apar).*(ps8.Color)+([1 1 1].*(1-apar))
ps0.Color=(apar).*(ps0.Color)+([1 1 1].*(1-apar))

set(gcf,'color','white')
grid on
xlabel('Disorder W(J)')
legend('U=0','U=2.7J','Data')
ylabel('\Delta x (sites)')
title('L=12 Comparison')

%{
figure(1)
subplot(1,2,2)
hold off
p0=plot(g2dU0Store(:,1),xisU0,'--')
hold on
ps8=plot(Thry12s(:,1),Thry12s(:,2),'-')
pe8=errorbar(Dat12s(:,1),Dat12s(:,2),Dat12s(:,3),'o')
pe8.Color=ps8.Color;
%}


