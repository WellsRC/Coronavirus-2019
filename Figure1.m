clear;
close all;
CC=[hex2rgb('#9BC01C');hex2rgb('#F5BE41');hex2rgb('#2D4262')];
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.1,0.64,0.283865546218487,0.341162790697676]); 
load('Daily_Prob_Expect.mat');
b=bar([minE:maxE],[MLExTS;MLExS-MLExTS;MLExTNS-(MLExS)]','stacked','LineStyle','none');

for ii=1:length(CC(:,1))
        b(ii).FaceColor = 'flat';
        b(ii).CData = CC(ii,:);
end
    
startDateofSim = datenum('12-06-2019');% Start date
XTL=datestr([startDateofSim+[0:4:(maxE-1)]],'mm-dd-yy');
yh=ylabel({'Expected number of infected','cases exported from China'},'Fontsize',18);
xlabel('Date','Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:4:maxE],'XTickLabel',XTL,'Xminortick','on','Yminortick','on');
xtickangle(45);
xlim([1 53.5])
ylim([0 8]);

legend({'Screening and travel ban','Screening and no travel ban','No screening and no travel ban'},'Location','NorthWest');
legend boxoff;
text(yh.Extent(1),max(ylim),'A','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.44,0.64,0.283865546218487,0.341162790697676]); 
plot([minE:maxE],MPTS,'color',CC(1,:),'LineWidth',2); hold on
plot([minE:maxE],MPS,'color',CC(2,:),'LineWidth',2); hold on
plot([minE:maxE],MPNS,'color',CC(3,:),'LineWidth',2);
yh=ylabel({'Probability infected','cases exported from China'},'Fontsize',18);
xlabel('Date','Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:4:maxE],'XTickLabel',XTL,'Xminortick','on','Yminortick','on');
xtickangle(45);
xlim([1 53])
ylim([0 1]);
text(yh.Extent(1),max(ylim),'B','Fontsize',32,'FontWeight','bold');



subplot('Position',[0.1,0.165,0.283865546218487,0.341162790697676]); 


plot([minE:maxE],MCPS,'color',CC(2,:),'LineWidth',2);
xlabel('Date of symptom onset','Fontsize',18);
yh=ylabel({'Probability at least one','one infected case exported','since 12-06-2019'},'Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:4:maxE],'XTickLabel',XTL,'Xminortick','on','Yminortick','on');
xtickangle(45);
xlim([1 53])
ylim([0 1]);
text(yh.Extent(1),max(ylim),'C','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.44,0.165,0.283865546218487,0.341162790697676]); 
[IncC,IncW,IncO]=IncidenceData;
I=cumsum(IncC(:,2)+IncW(:,2)+IncO(:,2));
load('Weighted_Travel_Inubation.mat');



contourf([0:52],ptravel.*w,MLE,[0:0.05:0.95 0.999],'LineStyle','none');
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'Xtick',[0:4:52],'XTicklabel',XTL,'Xminortick','on','Yminortick','on');
xtickangle(45);
xlabel('Date of symptom onset','Fontsize',18);
yhh=ylabel('Probability of travel per day','Fontsize',18);
y=colorbar;
yh=ylabel(y,{'Probability at least one','one infected case exported','since 12-06-2019'});
yh.Rotation=270;
yh.Position=[5.119364897410075,0.500000476837158,0];
y.Position=[0.730440414490399,0.165225744476465,0.011204481792717,0.341018245100029];
y.Limits=[0 1];
load('ColorM','CMC');
colormap(CMC)
text(yhh.Extent(1),max(ylim),'D','Fontsize',32,'FontWeight','bold');