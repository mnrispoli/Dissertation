%% axial lattice and 2D lattice heating

sc_2d=csvread('2DScatterWn.csv')';
sc_ax=csvread('AxialScatterWn.csv')';

V2D=[0.4 1 2 4 6 8]

ts=[0:2:12]
Vax=10.^(2*(0:0.1:0.4)).*(250/6.3096);

for vv=1:length(V2D)
    temp=find(V2D(vv)<sc_2d(:,1));
    Thry2d(vv)=2.*sc_2d(temp(1),2);
end

for vv=1:length(Vax)
    temp=find(Vax(vv)<sc_ax(:,1));
    Thryax(vv)=sc_ax(temp(1),2);
end




filenamesax={ 'TotalLossAxialV3.2.csv', ...
    'TotalLossAxialV3.3.csv', ...
    'TotalLossAxialV3.4.csv', ...
    'TotalLossAxialV3.5.csv', ...
    'TotalLossAxialV3.6.csv' }


for ff=1:length(filenamesax)
    data_ax(:,:,ff)=csvread(filenamesax{ff})
end


inds=1:length(ts)


           
%figure(1)
for ii=1:size(data_ax,3)
    %errorbar(data(:,1,ii),data(:,2,ii),data(:,3,ii),'o')
    %hold on
    
    f=fit(data_ax(inds,1,ii),data_ax(inds,2,ii),'exp1')
    fgs_ax(ii)=-f.b
    
    fgs_axErr(ii,:)=mean(abs((confint(f,0.6827))-ones(2,1)*coeffvalues(f)),1)
end


fitinds = 1:length(Vax)
x=Vax(fitinds);
y=fgs_ax(fitinds);
yax=y;
ybkg_ax = yax-Thryax;

fitfunc = @(b,x) 0.0305+b(2).*x.^b(3)

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,0],...
               'Upper',[1,100,1],...
               'StartPoint',[0.03 1 0.6]);
           
nlm=fitnlm(x,y,fitfunc,[0.0305 0.001 0.77])

cest = nlm.Coefficients.Estimate
cestErr = nlm.Coefficients.SE

vs=linspace(0,Vax(end),100);

figure(2)
plot(x,y,'o')
hold on
plot(vs,fitfunc(cest,vs),'--')




%%

filenames2d={ 'TotalLossAxialV3.2.csv', 
    'MiddleLoss2D1.0.csv', 
    'MiddleLoss2D2.0.csv', 
    'MiddleLoss2D4.0.csv',
    'MiddleLoss2D6.0.csv',
    'MiddleLoss2D8.0.csv'}

filenames2d={ 'TotalLossAxialV3.2.csv', 
    'TotalLoss2D1.0.csv', 
    'TotalLoss2D2.0.csv', 
    'TotalLoss2D4.0.csv',
    'TotalLoss2D6.0.csv',
    'TotalLoss2D8.0.csv'}

for ff=1:length(filenames2d)
    data(:,:,ff)=csvread(filenames2d{ff})
end

inds=4:length(ts)

%figure(1)
for ii=1:size(data,3)
    %errorbar(data(:,1,ii),data(:,2,ii),data(:,3,ii),'o')
    %hold on
    
    f=fit(data(inds,1,ii),data(inds,2,ii),'exp1')
    fgs_2d(ii)=-f.b
    
    fgs_2dErr(ii,:)=mean(abs((confint(f,0.6827))-ones(2,1)*coeffvalues(f)),1)
end

fitinds = 1:length(V2D)
x=V2D(fitinds);
y=fgs_2d(fitinds);
y2d=y;
ybkg_2d = y2d-Thry2d;
fitfunc = @(b,x) 0.0305+b(2).*x.^b(3)

nlm=fitnlm(x,y,fitfunc,[0.0315 1 0.6])

cest = nlm.Coefficients.Estimate
cestErr = nlm.Coefficients.SE

vs=linspace(0,12,100);

figure(3)
plot(x,y,'o')
hold on
plot(vs,fitfunc(cest,vs),'--')


%%
sc2dout=[V2D' y2d'-0.031 sqrt(fgs_2dErr(:,2).^2+0.003.^2)]
csvwrite('Sc2DOut.csv',sc2dout)
scaxout=[Vax' yax'-0.031 sqrt(fgs_axErr(:,2).^2+0.003.^2)]
csvwrite('ScAxOut.csv',scaxout)