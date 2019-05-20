clear all
ns=4;
basisExp=zeros(ns,ns);
xx=1:(ns-1);
xi=1;
ys=exp(-abs(xx)/xi);

basis=[1 0 1 0; 1 0 0 1; 0 1 1 0; 0 1 0 1;];

rvec=diag(rand(1,4).*2-1)*1;

uint=diag([ys(2) ys(3) ys(1) ys(2)]);

Hnon=diag(sum(basis*rvec,2));
Hint=diag(sum(uint,2));

yint=[1 1 1 1]./sqrt(4);

ts=linspace(0,10*2*pi,10000);
dt=ts(2)-ts(1);
yt_int(:,1)=yint';
yt_non(:,1)=yint';

for tt=2:length(ts)
    yt_int(:,tt)=expm(-i.*dt.*(Hnon+Hint))*yt_int(:,tt-1);
    yt_non(:,tt)=expm(-i.*dt.*(Hnon))*yt_non(:,tt-1);
    
    rhoa_int=(yt_int([1,2],tt)*yt_int([1,2],tt)')+ ...
        (yt_int([3,4],tt)*yt_int([3,4],tt)');
    [V,D]=eig(rhoa_int);
    Sint(tt)=-log(sum(abs(diag(D)).^2));
    
    rhoa_non=(yt_non([1,2],tt)*yt_non([1,2],tt)')+ ...
        (yt_non([3,4],tt)*yt_non([3,4],tt)');
    [V,D]=eig(rhoa_non);
    Snon(tt)=-log(sum(abs(diag(D)).^2));
end

figure(1)
plot(ts,Sint)
hold on
plot(ts,Snon)