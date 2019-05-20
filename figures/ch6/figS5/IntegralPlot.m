% for some of the ref replies we might generate some new plots


nf8sTh = csvread('nf_th_8s.csv')
nf12sTh = csvread('F_l1_12sites.csv');
nf8s = csvread('nf_8s_dat.csv')';
nf12s = csvread('nf_12s_dat.csv')';

inds = [1:8]
[ATZF12, TZF12, TZF12Err, DW12] = TrapzIntErr(nf12s)
[ATZF8, TZF8, TZF8Err, DW8] = TrapzIntErr(nf8s(inds,:))


indsDat = [1:8]
DFDat=[nf12s(:,1), nf12s(:,2)-nf8s(indsDat,2), sqrt(nf12s(:,3).^2+nf8s(indsDat,3).^2)]
[ADTZF, DTZF, DTZFErr, DDW] = TrapzIntErr(DFDat)

for kk=1:length(nf12sTh(:,1))
[ii,jj]=min(abs(nf8sTh(:,1)-nf12sTh(kk,1)));
inds(kk)=jj;
end


%{
[ATZF12_th, TZF12_th, TZF12Err_th, DW12_th]=TrapzIntErr(nf12sTh(3:end,:));
[ATZF8_th, TZF8_th, TZF8Err_th, DW8_th]=TrapzIntErr(nf8sTh(inds(3:end),:));
%}
for kk=1:length(nf8s(indsDat,1))
[ii,jj]=min(abs(nf8sTh(:,1)-nf8s(kk,1)));
inds8th(kk)=jj;
end
[ATZF8_thD, TZF8_thD, TZF8Err_thD, DW8_thD]=TrapzIntErr(nf8sTh(inds8th,:));

for kk=1:length(nf12s(:,1))
[ii,jj]=min(abs(nf12sTh(:,1)-nf12s(kk,1)));
inds12th(kk)=jj;
end
[ATZF12_thD, TZF12_thD, TZF12Err_thD, DW12_thD]=TrapzIntErr(nf12sTh(inds12th,:));
%{
[ATZF8_th_all, TZF8_th_all, TZF8Err_th_all, DW8_th_all]=TrapzIntErr(nf8sTh(:,:));
[ATZF12_th_all, TZF12_th_all, TZF12Err_th_all, DW12_th_all]=TrapzIntErr(nf12sTh(:,:));
%}
DFTh=[nf12sTh(inds12th,1), ...
    nf12sTh(inds12th,2)-nf8sTh(inds8th,2),  ...
    sqrt(0.*nf12sTh(inds12th,2).^2+0.*nf8sTh(inds8th,2).^2) ];

[ADTZF_th, DTZF_th, DTZFErr_th, DDW_th] = TrapzIntErr(DFTh);

figure(1)
subplot(1,2,1)
hold off

p8=plot(DW8_thD,TZF8_thD, ...
    'LineWidth',2);

hold on

p12=plot(DW12_thD,TZF12_thD, ...
        'LineWidth',2);

pe8=errorbar(DW8, TZF8, TZF8Err, 'o', ...
    'LineWidth',2, ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', 'White', ...
    'CapSize',0 ...
    );

pe12=errorbar(DW12, TZF12, TZF12Err, 's', ...
    'LineWidth',2, ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', 'White', ...
    'CapSize',0 ...
    );

pe12.Color=p12.Color
pe8.Color=p8.Color

apar=0.6;
p12.Color=(apar).*(p12.Color)+([1 1 1].*(1-apar))
p8.Color=(apar).*(p8.Color)+([1 1 1].*(1-apar))

grid on
xlabel('Depth W(J)')
ylabel('\int_{W_o}^{W} dw F_{L}(w)')
legend('L=8','L=12','location','Northwest')
xlim([3 10.5])

title('Integral of Fluctuations')

subplot(1,2,2)
hold off
pd = plot(DDW_th, DTZF_th)
hold on
pdErr = errorbar(DDW, DTZF, DTZFErr, 's')
pdErr.Color=pd.Color

hold on
plot([0 11],[0 0], 'k--')
grid on
xlabel('Depth W(J)')
ylabel('\int_{W_o}^{W}dw [F_{12}(w)-F_{8}(w)]')
title('Integral of Difference in Fluctuations')
xlim([3 10.5])



set(gcf,'color','white')






