close all;
load('Test_Travel_Ban.mat');

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.710882352941175,0.648680167762309,0.283865546218487,0.341162790697676]); 
histogram(UMLETS,[0:2:140],'Facecolor','b','LineStyle','none','FaceAlpha',0.3,'Normalization','probability'); hold on
histogram(UMLES,[0:2:140],'Facecolor','k','LineStyle','none','FaceAlpha',0.3,'Normalization','probability'); hold on
histogram(UMLE,[0:2:140],'Facecolor','r','LineStyle','none','FaceAlpha',0.3,'Normalization','probability'); hold on
xlabel('Expected number of cases','Fontsize',18);
yh=ylabel('Frequency','Fontsize',18);
legend('Travel restriction','No travel restriction','No Restrictions and no screening');
legend boxoff;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16);
box off;
xlim([0 140]);

startDateofSim = datenum('12-06-2019');% Start date

XTL=datestr([startDateofSim+[0:(length(IncO(:,2))-1)]],'mm-dd-yy');
subplot('Position',[0.065,0.648680167762309,0.588677688466691,0.341162790697676]); 
for ii=1:53
   patch(ii+[-0.35 -0.35 0.35 0.35],prctile(UMLExTS(:,ii),[2.5 97.5 97.5 2.5]),'b','Facealpha',0.3,'LineStyle','none'); hold on
   plot(ii+linspace(-0.35,0.35,101),MLExTS(ii).*ones(101,1),'b','LineWidth',2); hold on
   
   patch(ii+[-0.35 -0.35 0.35 0.35],prctile(UMLExS(:,ii),[2.5 97.5 97.5 2.5]),'k','Facealpha',0.3,'LineStyle','none'); hold on
   plot(ii+linspace(-0.35,0.35,101),MLExS(ii).*ones(101,1),'k','LineWidth',2); hold on
   
   patch(ii+[-0.35 -0.35 0.35 0.35],prctile(UMLEx(:,ii),[2.5 97.5 97.5 2.5]),'r','Facealpha',0.3,'LineStyle','none'); hold on
   plot(ii+linspace(-0.35,0.35,101),MLEx(ii).*ones(101,1),'r','LineWidth',2); hold on
end
xlim([0.5 length(IncO(:,2))+0.5]);
scatter([1:length(IncO(:,2))],IncO(:,2),40,'k','filled');
yh=ylabel({'Expected number of','cases missed'},'Fontsize',18);
xlabel('Date of symptom onset','Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:length(IncO(:,2))],'XTickLabel',XTL,'YTick',[0:5:30],'Yminortick','on');
xtickangle(45);
ylim([0 15]);

subplot('Position',[0.065,0.648680167762309,0.588677688466691,0.341162790697676]); 