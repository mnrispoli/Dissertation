rnStore=[];
rnmStore=[];

for rep=1:1000
rvec=rand([1,1000]);
rsort=sort(rvec);
gaps=rsort(2:end)-rsort(1:end-1);

for gg=1:length(gaps)-1
rn(gg)=min(gaps(gg),gaps(gg+1))./max(gaps(gg),gaps(gg+1));
end

rmat=rand(1000);
rsort=eig(rmat+ctranspose(rmat));
gaps=rsort(2:end)-rsort(1:end-1);

for gg=1:length(gaps)-1
rnm(gg)=min(gaps(gg),gaps(gg+1))./max(gaps(gg),gaps(gg+1));
end

rnStore=[rnStore rn];
rnmStore=[rnmStore rnm];

end
%%
figure(1)
histogram(rnStore,50,'Normalization','probability','EdgeColor',[42 157 242]/256,'FaceColor',[103 175 229]/256);
hold on
histogram(rnmStore,50,'Normalization','probability','EdgeColor',[252 170 71]/256,'FaceColor',[252 183 101]/256);
set(gcf,'color','white')
xlim([0 1])
ylim([0 0.045])
xlabel('r_n')
ylabel('P(r_n)')
%%

[rnvals,rnbins]=histcounts(rnStore,50,'Normalization','probability')
[rnmvals,rnmbins]=histcounts(rnmStore,50,'Normalization','probability')

csvwrite('poissonstats.csv',[rnbins(2:end)-0.01; rnvals]')
csvwrite('goestats.csv',[rnmbins(2:end)-0.01; rnmvals]')