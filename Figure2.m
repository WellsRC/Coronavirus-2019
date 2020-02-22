%% Country level risk assesment
close all;
load('Probability_Travel_Infection_6733.mat','F','pc');
ptravel=pc(F==max(F));
figure('units','normalized','outerposition',[0 0 1 1]);

[IncC,IncW,IncO]=IncidenceData;
I=cumsum(IncC(:,2)+IncW(:,2)+IncO(:,2));
load('Weighted_Travel_Infectious_6733.mat');
minE=-22;
maxE=72;

startDateofSim = datenum('12-06-2019');% Start date
XTL=datestr([startDateofSim+[0:1:(maxE-1)]],'mm-dd-yy');

INDX=datenum('01-23-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
INDX2=datenum('01-25-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero

contourf([(minE+1):maxE],ptravel.*w,log10(abs(MLE(:,2:end)-MLE(:,1:end-1))),[-5:1:-1],'LineStyle','none');
box off;
xlim([1 maxE]);
ylim([0 0.5*10^(-3)]);
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'Xtick',[1:1:maxE],'XTicklabel',XTL,'Xminortick','on','Yminortick','on','YTick',[0:0.2:2]*10^(-4));
xtickangle(45);
xlabel('Date','Fontsize',18);
%title({'Probability of exportation since start of outbreak'},'Fontsize',18);%
yhh=ylabel('Country level probability of travel per day','Fontsize',18);
y=colorbar;
yhC=ylabel(y,{'Risk of initial importation event','(Probability)'});
yhC.Rotation=270;
yhC.Position=[5.88,mean([[-5 -1]]),0];
y.Position=[0.913738733818127,0.133738601823708,0.011204481792717,0.79128672745694];
y.Limits=[-5 -1];
y.Ticks=[-5:-1];
y.TickLabels={'10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}'};
load('ColorM','CMC');
colormap(CMC)


hold on
xlim([30 57])
ylim([0 2*10^(-4)])
load('Weight_Flights');
Cor=zeros(9,2);
for ii=1:length(FC)
    if(FC{ii,2}>=INDX)
        tf = strcmp({FC{ii,1}},{Flight_NW{:,1}});
        scatter(FC{ii,2}-1,Flight_NW{tf,2}*ptravel,25,'k','filled'); hold on;
        text((FC{ii,2}-1)+0.06,(Flight_NW{tf,2}*ptravel)+0.0000015,15,{FC{ii,1}},'Fontsize',14,'Rotation',45);
        Cor(ii,:)=[Flight_NW{tf,2} (FC{ii,2}-1)];
    else
        tf = strcmp({FC{ii,1}},{FlightAll{:,1}});
        scatter(FC{ii,2}-1,FlightAll{tf,2}*ptravel,25,'k','filled'); hold on;
        text((FC{ii,2}-1)+0.06,(FlightAll{tf,2}*ptravel)+0.0000015,15,{FC{ii,1}},'Fontsize',14,'Rotation',45);
        Cor(ii,:)=[FlightAll{tf,2} (FC{ii,2}-1)];
    end
end
pp=plot([INDX INDX],[0 2*10^(-4)],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
plot([INDX2 INDX2],[0 2*10^(-4)],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
legend(pp,'Enforcement of lockdown');
legend boxoff;
[r,p]=corr(Cor)
round(p(2,1),3)
%text(30.5,4.4*10^(-4)*0.95,['r=' num2str(round(r(1,2),2)) ' (p=' num2str(round(p(1,2),3)) ')'],'Fontsize',16);
clear;