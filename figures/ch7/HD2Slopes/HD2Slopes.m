folderprefix='C:\Users\Matthew\Documents\RbPapers\MBL_Critical_Behavior\data_repo\';
foldername={'W2_time_trace', ...
            'W8.4_time_trace', ...
            'W9.6_time_trace', ...
            'W10.3_time_trace', ...
            'W11.5_time_trace', ...
            'W13.6_time_trace', ...
            'W15.6_time_trace', ...
            'W17.8_time_trace'};
filename = '\HD2Slopes_data.csv';
HDfilename = '\HD2_data.csv';
HDThfilename = '\HD2_theory.csv';

DataSlopes=zeros(length(foldername),3);

for ff=1:length(foldername)
    A=csvread([folderprefix foldername{ff} filename]);
    DataSlopes(ff,:)=abs(A);
    D{ff}=csvread([folderprefix foldername{ff} HDfilename]);
    temp=D{ff};
    csvwrite(['HD2_data_' foldername{ff} '.csv'],temp(3:end,:));
    Dth{ff}=csvread([folderprefix foldername{ff} HDThfilename]);
    temp=Dth{ff};
    csvwrite(['HD2_thry_' foldername{ff} '.csv'],temp);
end

csvwrite('HD2Slopes_data.csv',DataSlopes);
%HD2Slopes_data

B=csvread('C:\Users\Matthew\Documents\RbPapers\MBL_Critical_Behavior\data_repo\W17.8_time_trace\HD2Slopes_Theory.csv');
csvwrite('HD2Slopes_Theory.csv',abs(B))
%csvwrite('HDSlopes_Theory.csv',abs(B))

figure(33)
plot(B(:,1),abs(B(:,2)))
hold on
errorbar(DataSlopes(:,1),abs(DataSlopes(:,2)),DataSlopes(:,3),'o')
ylim([0 4])
xlim([0 25])