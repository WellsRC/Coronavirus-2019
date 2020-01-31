close all;
load('Test_Travel_Ban.mat');
C1=sum(UxT,2);
C2=sum(UxTN,2);
C3=sum(UxTNS,2);
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.710882352941175,0.648680167762309,0.283865546218487,0.341162790697676]); 
histogram(C1,[0:5:450],'Facecolor','b','LineStyle','none','FaceAlpha',0.3,'Normalization','probability'); hold on
histogram(C2,[0:5:450],'Facecolor','k','LineStyle','none','FaceAlpha',0.3,'Normalization','probability'); hold on
histogram(C3,[0:5:450],'Facecolor','r','LineStyle','none','FaceAlpha',0.3,'Normalization','probability'); hold on
xlabel('Expected number of cases','Fontsize',18);
yh=ylabel('Frequency','Fontsize',18);
legend('Travel restriction','No travel restriction','No Restrictions and no screening');
legend boxoff;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16);
box off;
xlim([0 450]);

startDateofSim = datenum('12-29-2020');% Start date

XTL=datestr([startDateofSim+[0:(length(Inc(:,2))-1)]],'mm-dd-yy');
subplot('Position',[0.065,0.648680167762309,0.588677688466691,0.341162790697676]); 
for ii=1:30
   patch(ii+[-0.35 -0.35 0.35 0.35],prctile(UxT(:,ii),[2.5 97.5 97.5 2.5]),'b','Facealpha',0.3,'LineStyle','none'); hold on
   plot(ii+linspace(-0.35,0.35,101),MLExT(ii).*ones(101,1),'b','LineWidth',2); hold on
   
   patch(ii+[-0.35 -0.35 0.35 0.35],prctile(UxTN(:,ii),[2.5 97.5 97.5 2.5]),'k','Facealpha',0.3,'LineStyle','none'); hold on
   plot(ii+linspace(-0.35,0.35,101),MLExTN(ii).*ones(101,1),'k','LineWidth',2); hold on
   
   patch(ii+[-0.35 -0.35 0.35 0.35],prctile(UxTNS(:,ii),[2.5 97.5 97.5 2.5]),'r','Facealpha',0.3,'LineStyle','none'); hold on
   plot(ii+linspace(-0.35,0.35,101),MLExTNS(ii).*ones(101,1),'r','LineWidth',2); hold on
end
xlim([0.5 length(Inc(:,2))+0.5]);
scatter([1:length(IncO(:,2))],IncO(:,2),40,'k','filled');
yh=ylabel({'Expected number of','cases missed'},'Fontsize',18);
xlabel('Date of symptom onset','Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:length(Inc(:,2))],'XTickLabel',XTL,'YTick',[0:10:100],'Yminortick','on');
xtickangle(45);
ylim([0 100]);