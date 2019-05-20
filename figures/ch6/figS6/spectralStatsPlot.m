clear all
close all
clc

prefix=''
RstatL5=csvread([prefix 'RStat_L5.csv'])
RstatL6=csvread([prefix 'RStat_L6.csv'])
RstatL7=csvread([prefix 'RStat_L7.csv'])
RstatL8=csvread([prefix 'RStat_L8.csv'])




figure(1)

hold off
ps6=plot(RstatL6(:,1)./2,RstatL6(:,2),'-', ...
    'LineWidth',2)
hold on
ps8=plot(RstatL8(:,1)./2,RstatL8(:,2),'-', ...
    'LineWidth',2)




set(gcf,'color','white')
grid on
xlabel('Disorder W(J)')
legend('L=6','L=8')
ylabel('\langle r \rangle')
title('Spectral Statistics')

figure(2)

hold off
ps6=plot(RstatL7(:,1)./2,RstatL7(:,2),'-', ...
    'LineWidth',2)
hold on
ps8=plot(RstatL8(:,1)./2,RstatL8(:,2),'-', ...
    'LineWidth',2)




set(gcf,'color','white')
grid on
xlabel('Disorder W(J)')
legend('L=7','L=8')
ylabel('\langle r \rangle')
title('Spectral Statistics')

%{
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
ps6=plot(RstatL5(:,1),RstatL5(:,2),'k--', ...
    'LineWidth',2)
hold on
ps8=plot(RstatL7(:,1),RstatL7(:,2),'-', ...
    'LineWidth',2)


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
    %}