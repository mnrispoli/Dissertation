clear all
close all
Ns=31;
x=linspace(0,Ns,1001);
x=x(1:end-1);
dx=x(2)-x(1);
hv=10.*diag(cos(pi.*x).^2);
hpv=(zeros(size(x)));
hpv(:,1:2)=[2 -1];
hp=toeplitz(hpv);
hp(1,end)=-1;
hp(end,1)=-1;

H=hp./(dx.^2)+hv;
[V,D]=eig(H);
%%
figure(1)
%subplot(1,2,1)
%plot(diag(D));

%subplot(1,2,2)
hold off
plot(diag(hv)./max(diag(hv)))
hold on
for ee=1:Ns
    if mod(ee,2)==1
        pree(:,ee)=((V(:,ee)))./sign(V(500,ee));
    else
        pree(:,ee)=i.*V(:,ee)./sign(V(505,ee));
    end
    %pree(:,ee)=pree(:,ee)./max(abs(pree(:,ee)));
end
%prw=abs(sum(V,2)).^2;
%prw=prw./max(prw);
%plot(abs(pree(:,1:3)))
hold on
plot(abs(sum(pree(:,:),2)).^2)
%plot(prw)

%%
%{
pq=zeros(size(V(:,1)))
for qq=1:100
    pq=pq+V(:,1).*(exp(-i.*pi.*[-500:500].*(qq-1)./100))';
end
%}