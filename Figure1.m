clear;
clc;
CC=[hex2rgb('#2D4262');hex2rgb('#F5BE41');];

CCT=[hex2rgb('#20948B')];


INDX=datenum('01-23-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
INDX2=datenum('01-25-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
INDXMV=datenum('01-1-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero

load('Daily_Prob_Expect_6733.mat');

startDateofSim = datenum('12-06-2019');% Start date
XTL=datestr([startDateofSim+[0:4:(maxE-1)]],'mm-dd-yy');

figure('units','normalized','outerposition',[0 0 1 1]);

subplot('Position',[0.045,0.165,0.283865546218487,0.3]); 

b=bar([minE:maxE],[MLExS;MLExNS-(MLExS)]','stacked','LineStyle','none'); hold on;

for ii=1:2
        b(ii).FaceColor = 'flat';
        b(ii).CData = CC(ii,:);
end
    
xlabel('Date','Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:4:maxE],'XTickLabel',XTL,'Xminortick','on','Yminortick','on');

xtickangle(45);
xlim([1 maxE+0.5])
ylim([0 35]);
title('No travel lockdown','Fontsize',18);
legend(b,{'Incubation','Symptomatic'},'Location','NorthWest');
legend boxoff;
yh=ylabel({'Number of exported cases'},'Fontsize',18);

text(yh.Extent(1),max(ylim)*1.1,'D','Fontsize',32,'FontWeight','bold');


subplot('Position',[0.367+0.01,0.165,0.283865546218487,0.3]); 

b=bar([minE:maxE],[MLExTS;MLExTNS-(MLExTS)]','stacked','LineStyle','none'); hold on

for ii=1:2
        b(ii).FaceColor = 'flat';
        b(ii).CData = CC(ii,:);
end
    
xlabel('Date','Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:4:maxE],'XTickLabel',XTL,'Xminortick','on','Yminortick','on');
plot([INDX INDX],[0 35],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
plot([INDX2 INDX2],[0 35],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
xtickangle(45);
xlim([1 maxE+0.5])
ylim([0 35]);
title('Travel lockdown','Fontsize',18);
legend(b,{'Incubation','Symptomatic'},'Location','NorthWest');
legend boxoff;
yh=ylabel({'Number of exported cases'},'Fontsize',18);

text(yh.Extent(1),max(ylim)*1.1,'E','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.689+0.02,0.165,0.283865546218487,0.3]); 

b=bar([minE:maxE],[MLExS-MLExTS;(MLExNS-(MLExTNS))-(MLExS-MLExTS)]','stacked','LineStyle','none'); hold on;

for ii=1:2
        b(ii).FaceColor = 'flat';
        b(ii).CData = CC(ii,:);
end
    
xlabel('Date','Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:4:maxE],'XTickLabel',XTL,'Xminortick','on','Yminortick','on');
plot([INDX INDX],[0 35],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
plot([INDX2 INDX2],[0 35],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
xtickangle(45);
xlim([1 maxE+0.5])
ylim([0 35]);
title('Cases averted by travel lockdown','Fontsize',18);
legend(b,{'Incubation','Symptomatic'},'Location','NorthWest');
legend boxoff;
yh=ylabel({'Number of exported cases'},'Fontsize',18);

text(yh.Extent(1),max(ylim)*1.1,'F','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.045,0.64,0.283865546218487,0.3]); 
%plot([minE:(INDX)],MPTS(1:(1+INDX-minE)),'color',CC(1,:),'LineWidth',2); hold on
plot([minE:(INDX)],MPTNS(1:(1+INDX-minE)),'color',CCT(1,:),'LineWidth',2); hold on


%plot([(INDX):maxE],MPTS((1+INDX-minE):end),'-.','color',CC(1,:),'LineWidth',2); hold on
plot([(INDX):maxE],MPTNS((1+INDX-minE):end),'-.','color',CCT(1,:),'LineWidth',2); 

%\\% LB=prctile(UMPTNS,2.5);
%\\% UB=flip(prctile(UMPTNS,97.5));
%\\% patch([[minE:maxE] flip([minE:maxE])],[LB UB],CCT(1,:),'LineStyle','none','Facealpha',0.35);

yhB=ylabel({'Probability'},'Fontsize',18);
xlabel('Date','Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:4:maxE],'XTickLabel',XTL,'Xminortick','on','Yminortick','on');
plot([INDX INDX],[0 1],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
plot([INDX2 INDX2],[0 1],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
xtickangle(45);
xlim([1 maxE])
ylim([0 1]);
title('Daily risk of exportation','Fontsize',18);
text(yhB.Extent(1),max(ylim)*1.1,'A','Fontsize',32,'FontWeight','bold');

%legend({'Incubation','Infected'},'Location','NorthWest');
%legend boxoff;

subplot('Position',[0.367+0.01,0.64,0.283865546218487,0.3]); 

%plot([minE:maxE],MCPTS,'color',CC(1,:),'LineWidth',2); hold on
plot([minE:(INDX)],MCPTNS(1:(1+INDX-minE)),'color',CCT(1,:),'LineWidth',2); hold on
plot([(INDX):maxE],MCPTNS((1+INDX-minE):end),'-.','color',CCT(1,:),'LineWidth',2); hold on

%\\% LB=prctile(UMCPTNS,2.5);
%\\% UB=flip(prctile(UMCPTNS,97.5));
%\\% patch([[minE:maxE] flip([minE:maxE])],[LB UB],CCT(1,:),'LineStyle','none','Facealpha',0.35);
plot([INDX INDX],[0 1],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
plot([INDX2 INDX2],[0 1],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);

yhB=ylabel({'Probability'},'Fontsize',18);
xlabel('Date','Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:4:maxE],'XTickLabel',XTL,'Xminortick','on','Yminortick','on');
xtickangle(45);
xlim([1 maxE])
ylim([0 1]);
title({'Cumulative risk of exportation'},'Fontsize',18);
text(yhB.Extent(1),max(ylim)*1.1,'B ','Fontsize',32,'FontWeight','bold');

% legend({'Incubation','Infected'},'Location','NorthWest');
% legend boxoff;

subplot('Position',[0.689+0.02,0.64,0.283865546218487,0.3]); 
TT=MCPTNS(2:end)-MCPTNS(1:end-1);

plot([(minE+1):(INDX)],TT(1:(1+INDX-(minE+1))),'color',CCT(1,:),'LineWidth',2); hold on
plot([(INDX):maxE],TT((1+INDX-(minE+1)):end),'-.','color',CCT(1,:),'LineWidth',2); hold on

plot([INDX INDX],[0 0.1],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
plot([INDX2 INDX2],[0 0.1],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);

%\\% LB=prctile(UMCPTNS(:,2:end)-UMCPTNS(:,1:end-1),2.5);
%\\% UB=flip(prctile(UMCPTNS(:,2:end)-UMCPTNS(:,1:end-1),97.5));
%\\% patch([[(minE+1):maxE] flip([(minE+1):maxE])],[LB UB],CCT(1,:),'LineStyle','none','Facealpha',0.35);

yhB=ylabel({'Probability'},'Fontsize',18);
xlabel('Date','Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:4:maxE],'XTickLabel',XTL,'Xminortick','on','Yminortick','on','YTick',[0:0.01:0.1]);
xtickangle(45);
xlim([1 maxE])
ylim([0 0.1]);
title({'Risk of initial exportation event'},'Fontsize',18);
text(yhB.Extent(1),max(ylim)*1.1,'C','Fontsize',32,'FontWeight','bold');
clear;