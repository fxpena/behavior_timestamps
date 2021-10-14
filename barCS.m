function barCS(x,ydata)
%Plots bars with error bars
%INPUTS
%x: scalar
%ydata: matrix where column 1 is CS+, column 2 is CS-
%   optionally, column 3 can be baseline
%NOTE: assumes a figure is already open
%for each column calculate the mean and SEM

%if there's more than 1 row of data (more than 1 mouse/test)
if size(ydata,1)>1
    m=mean(ydata);
    SEM=std(ydata)/sqrt(size(ydata,1)); %use the # of rows
else %for an indivdual mouse
    m=ydata;
    SEM=zeros(size(m));
end
%first place the bars
b=bar(x,m,'FaceColor',[0 0.7 0.6],'BarWidth',0.9,'LineWidth',1.5); %teal bars
%change the second series to yellow
b(2).FaceColor=[0.9 0.9 0]; %the second series is the CS-
switch size(ydata,2)
    case 2 %if there are 2 columns
        %determine coordinates for the errorbars
        x2=[x-0.15 x+0.15]; %each column slightly offset from the center
        %set LineSpec to black, crosses, no line between the means
        errorbar(x2,m,SEM,'kx','LineWidth',1.5)
    case 3 %if there are 3 columns
        %change the third series to white bars
        b(3).FaceColor=[1 1 1];
        %determine coordinates for the errorbars
        x2=[x-0.22 x x+0.22]; %each column slightly offset from the center
        %set LineSpec to black, crosses, no line between the means
        errorbar(x2,m,SEM,'kx','LineWidth',1.5)
    case 5 %if there are 5 columns
        %change colors of the other series
        b(4).FaceColor=[1 1 1]; %white
        b(3).FaceColor='magenta';
        b(5).FaceColor='blue';
        %determine coordinates for the errorbars
        x2=[x-0.31 x-0.155 x x+0.155 x+0.31]; %each column slightly offset from the center
        %set LineSpec to black, crosses, no line between the means
        errorbar(x2,m,SEM,'kx','LineWidth',1.5)
    case 7 %if there are 7 columns
        %change colors of the other series
        b(4).FaceColor=[1 1 1]; %white
        b(3).FaceColor='magenta';
        b(5).FaceColor='green';
        b(6).FaceColor='yellow';
        b(7).FaceColor='blue';
        %determine coordinates for the errorbars
        x2=[x-0.34 x-0.23 x-0.11 x x+0.11 x+0.23 x+0.34];
        %set LineSpec to black, crosses, no line between the means
        errorbar(x2,m,SEM,'kx','LineWidth',1.5)
end
end