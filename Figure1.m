clear;
clc;
close all;
CC=[hex2rgb('#9BC01C');hex2rgb('#F5BE41');hex2rgb('#2D4262')];




load('Daily_Prob_Expect.mat');

startDateofSim = datenum('12-06-2019');% Start date
XTL=datestr([startDateofSim+[0:4:(maxE-1)]],'mm-dd-yy');

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.1,0.165,0.283865546218487,0.341162790697676]); 


plot([minE:maxE],MCPS,'color',CC(2,:),'LineWidth',2);
xlabel('Date','Fontsize',18);
yhC=ylabel({'Probability at least one','one infected case exported','since start of outbreak'},'Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:4:maxE],'XTickLabel',XTL,'Xminortick','on','Yminortick','on');
xtickangle(45);
xlim([1 maxE])
ylim([0 1]);
text(yhC.Extent(1),max(ylim),'C','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.1,0.64,0.283865546218487,0.341162790697676]); 

b=bar([minE:maxE],[MLExTS;MLExS-MLExTS;MLExTNS-(MLExS)]','stacked','LineStyle','none');

for ii=1:length(CC(:,1))
        b(ii).FaceColor = 'flat';
        b(ii).CData = CC(ii,:);
end
    
ylabel({'Expected number of infected','cases exported from China'},'Fontsize',18);
xlabel('Date','Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:4:maxE],'XTickLabel',XTL,'Xminortick','on','Yminortick','on');
xtickangle(45);
xlim([1 maxE+0.5])
ylim([0 62]);

legend({'Screening and travel ban','Screening and no travel ban','No screening and no travel ban'},'Location','NorthWest');
legend boxoff;
text(yhC.Extent(1),max(ylim),'A','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.44,0.64,0.283865546218487,0.341162790697676]); 
plot([minE:maxE],MPTS,'color',CC(1,:),'LineWidth',2); hold on
plot([minE:maxE],MPS,'color',CC(2,:),'LineWidth',2); hold on
plot([minE:maxE],MPNS,'color',CC(3,:),'LineWidth',2);
yhB=ylabel({'Probability infected','cases exported from China'},'Fontsize',18);
xlabel('Date','Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:4:maxE],'XTickLabel',XTL,'Xminortick','on','Yminortick','on');
xtickangle(45);
xlim([1 maxE])
ylim([0 1]);
text(yhB.Extent(1),max(ylim),'B','Fontsize',32,'FontWeight','bold');

fprintf(['Expected date of exportation:' datestr([startDateofSim+([minE:maxE]*MPTS'./sum(MPTS)-1)],'mmm-dd-yy') '\n'])

subplot('Position',[0.44,0.165,0.283865546218487,0.341162790697676]); 
[IncC,IncW,IncO]=IncidenceData;
I=cumsum(IncC(:,2)+IncW(:,2)+IncO(:,2));
load('Weighted_Travel_Inubation.mat');



contourf([minE:maxE],ptravel.*w,MLE,[0:0.05:0.95 0.999],'LineStyle','none');
box off;
xlim([1 maxE]);
ylim([0 0.5*10^(-3)]);
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'Xtick',[1:4:53],'XTicklabel',XTL,'Xminortick','on','Yminortick','on','YTick',[0:5]*10^(-4));
xtickangle(45);
xlabel('Date','Fontsize',18);
yhh=ylabel('Probability of travel per day','Fontsize',18);
y=colorbar;
yhC=ylabel(y,{'Probability at least one','one infected case exported','since start of outbreak'});
yhC.Rotation=270;
yhC.Position=[6.262221881321509,0.500000476837158,0];
y.Position=[0.730440414490399,0.165225744476465,0.011204481792717,0.341018245100029];
y.Limits=[0 1];
load('ColorM','CMC');
colormap(CMC)
text(yhB.Extent(1),max(ylim),'D','Fontsize',32,'FontWeight','bold');

hold on

load('Weight_Flights');
Cor=zeros(9,2);
for ii=1:length(FC)
    tf = strcmp({FC{ii,1}},{FlightAll{:,1}});
    scatter(FC{ii,2}-1,FlightAll{tf,2}*ptravel,15,'k','filled'); hold on;
    if((ii==6)||(ii==5))
        text((FC{ii,2}-1)-0.5,(FlightAll{tf,2}*ptravel),15,{FC{ii,1}},'Fontsize',11,'HorizontalAlignment','right');
    else
        text((FC{ii,2}-1)+0.05,(FlightAll{tf,2}*ptravel)+0.000005,15,{FC{ii,1}},'Fontsize',11,'Rotation',45);
    end
    Cor(ii,:)=[FlightAll{tf,2} (FC{ii,2}-1)];
end
[r,p]=corr(Cor);
text(1.5,max(ylim)*0.95,['r=' num2str(round(r(1,2),2)) ' (p=' num2str(round(p(1,2),3)) ')'],'Fontsize',14);