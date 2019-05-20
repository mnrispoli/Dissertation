


% modified colorscale 04/15
c = [1 1 1; 0.12 0.25 1; 1/1.5 0.125/1.5 1/1.5; 1 0 0];
x = [0 0.25 0.75 1];
xx = linspace(0,1,100);
cc1 = interp1(x, c(:,1), xx);
cc2 = interp1(x, c(:,2), xx);
cc3 = interp1(x, c(:,3), xx);
alexmap = [cc1' cc2' cc3'];
pmax = [0.06 0.11 0.23 1.1];
ylimit = [0.28 0.28 0.4 1.3];


figure(2)
subplot(1,2,1)
J=0.011
tj=allUnique.*2.*pi.*J;
si=[-10:1:10];

imagesc(si,tj,yy(2:end,:))
axis('square')
caxis([0 1])
colormap(alexmap)
set(gcf,'color','white')
colorbar()
title('Data')
xlabel('Sites')
ylabel('time (\tau)')

Jeff=2*pi*fitted.J/1000;
%Feff=fitted2.F;

for tji=2:length(tj)
yeff(tji-1,:)=abs(besselj(abs(si), 4*pi*Jeff*allUnique(tji))).^2;
end

subplot(1,2,2)
imagesc(si,tj,yeff)
axis('square')
colormap(alexmap)
caxis([0 1])
colorbar()
title('Theory')
xlabel('Sites')
ylabel('time (\tau)')