%% Pan Equation Modeling
% Models the concentration of cocaine in the brain from infusions.
% 
% Inputs: 
% 	excel file or matrix table of infusions size(N animals, T timepoints)
% Outputs:
% 	corresponding model tracings of response to infusions as modeled by the Pan Equaition
%
% Author: Brian Carson and Sami Kummer
% Change Log
% 3/2018: Created by Brian Carson for mouse cocaine self administration data.
clear all, clc;
filename='ITIs for matlab 7-11-22.xlsx';
table=readmatrix(filename);
firstRow=1; % first row of data
firstColumn=1; % first coulmn of data
table=table(firstRow:end,firstColumn:end);
%%
d=.5;
a=9.637; %60;
alph=.642/60;
beta=.057/60;

base=0;
figure
hold on
col='rgb'; % colors for plotting each group
grp=[1 1 1 2 2 2 2 2 3 3 3]; % group definition
for j=1:size(table,1)
    ts1=find(diff(table(j,:)));%./60;


    x=zeros(size(table,2),1);
    for i=1:length(ts1)
        if ts1(i)>30
           % base=30;
        end
        c=[];
        t=ts1(i)-base:size(table,2);%/60;
        t=t-ts1(i);
        c=d*a*(exp(-beta*t)-exp(-alph*t));
        x(round((ts1(i)))-base:end)=x(round(ts1(i))-base:end)+c';
    end
    y(j,:)=x; % save this for plotting average of each group
    plot(1:7200,smooth(x,1),col(grp(j)))
    plot(ts1,ones(length(ts1),1)*60+2*j,'.','Color',col(grp(j)))
end
plot(mean(y(1:3,:)),col(1),'LineWidth',5) % select rows of each group for averages
plot(mean(y(4:8,:)),col(2),'LineWidth',5)
plot(mean(y(9:11,:)),col(3),'LineWidth',5)

