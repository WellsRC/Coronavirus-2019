%% Country level risk assesment
close all;
load('Probability_Travel_Infection.mat','F','pc');
ptravel=pc(F==max(F));
figure('units','normalized','outerposition',[0 0 1 1]);

[IncC,IncW,IncO]=IncidenceData;
I=cumsum(IncC(:,2)+IncW(:,2)+IncO(:,2));
load('Weighted_Travel_Infectious.mat');
minE=-22;
maxE=57;

startDateofSim = datenum('12-06-2019');% Start date
XTL=datestr([startDateofSim+[0:1:(maxE-1)]],'mm-dd-yy');

contourf([(minE+1):maxE],ptravel.*w,MLE(:,2:end)-MLE(:,1:end-1),[0:0.001:0.085],'LineStyle','none');
box off;
xlim([1 maxE]);
ylim([0 0.5*10^(-3)]);
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'Xtick',[1:1:maxE],'XTicklabel',XTL,'Xminortick','on','Yminortick','on','YTick',[0:7]*10^(-5));
xtickangle(45);
xlabel('Date','Fontsize',18);
%title({'Probability of exportation since start of outbreak'},'Fontsize',18);%
yhh=ylabel('Probability of travel per day','Fontsize',18);
y=colorbar;
yhC=ylabel(y,{'Likelihood of first exportation event'});
yhC.Rotation=270;
yhC.Position=[4.452698026384615,mean([0 0.085]),0];
y.Position=[0.913738733818127,0.133738601823708,0.011204481792717,0.79128672745694];
y.Limits=[0 0.085];
y.Ticks=[0:0.01:0.08];
load('ColorM','CMC');
colormap(CMC)


hold on
xlim([27 57])
ylim([0 7.5*10^(-5)])
load('Weight_Flights');
Cor=zeros(9,2);
for ii=1:length(FC)
    tf = strcmp({FC{ii,1}},{FlightAll{:,1}});
    scatter(FC{ii,2}-1,FlightAll{tf,2}*ptravel,15,'k','filled'); hold on;
%      if((ii==20))
%          text((FC{ii,2}-1)-0.2,(FlightAll{tf,2}*ptravel),15,{FC{ii,1}},'Fontsize',11,'HorizontalAlignment','right');
% %     elseif(ii==9)        
% %         text((FC{ii,2}-1)+0.015,(FlightAll{tf,2}*ptravel),15,{FC{ii,1}},'Fontsize',11);
%      else
        text((FC{ii,2}-1)+0.015,(FlightAll{tf,2}*ptravel)+0.00000075,15,{FC{ii,1}},'Fontsize',14,'Rotation',45);
%      end
    Cor(ii,:)=[FlightAll{tf,2} (FC{ii,2}-1)];
end
[r,p]=corr(Cor);
text(27.5,7.5*10^(-5)*0.95,['r=' num2str(round(r(1,2),2)) ' (p=' num2str(round(p(1,2),3)) ')'],'Fontsize',16);
clear;