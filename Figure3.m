close all;
clear;
CIP=hex2rgb('#2D4262');
CC=[hex2rgb('#F5BE41');hex2rgb('#CB0000');];
figure('units','normalized','outerposition',[0 0 1 1]);
load('Time_Screening.mat');
subplot('Position',[0.367+0.01,0.653,2.*0.283865546218487,0.3]); 
plot([1:14],MLET,'color',CIP(1,:),'LineWidth',2); hold on
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[0:14],'Xminortick','on','Yminortick','on');

LB=prctile(UMLET,2.5);
UB=flip(prctile(UMLET,97.5));
patch([[1:14] flip([1:14])],[LB UB],CIP,'LineStyle','none','Facealpha',0.35);
box off
xlim([1 14])
ylim([0 1]);
yh=ylabel({'Probability of identification'},'Fontsize',18);
xlabel('Days since exposure (infection)','Fontsize',18);

text(-0.66,max(ylim)*1.1,'B','Fontsize',32,'FontWeight','bold');

load('TravelDuringInfection.mat')


subplot('Position',[0.367+0.01,0.21,2.*0.283865546218487,0.3]); 
plot([0:14],(MLE-MLECT)./MLE,'color',CIP(1,:),'LineWidth',2); hold on;

LB=prctile((L-LCT)./L,2.5);
UB=flip(prctile((L-LCT)./L,97.5));
patch([[0:14] flip([0:14])],[LB UB],CIP,'LineStyle','none','Facealpha',0.35);
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[0:14],'Xminortick','on','Yminortick','on');
box off
xlim([0 14])
ylim([0 1]);
yh=ylabel({'Reduction in the probability of',' travel during incubation period'},'Fontsize',18);
xlabel('Days from infection to quarantine','Fontsize',18);

text(yh.Extent(1),max(ylim)*1.1,'C','Fontsize',32,'FontWeight','bold');


figure('units','normalized','outerposition',[0 0 1 1]);

subplot('Position',[0.367+0.01,0.21,2.*0.283865546218487,0.3]); 
load('Time_After_Arrival.mat','UMLE','UMLET','MLE','MLET')

% MLET is the time from arrival to symptom onset
histogram(UMLET,100,'LineStyle','none','FaceAlpha',0.3,'FaceColor',CC(1,:)); hold on
histogram(UMLE+MLET,100,'LineStyle','none','FaceAlpha',0.3,'FaceColor',CC(2,:)); hold on
box off;
ax1=gca;
ax1.YAxis.Visible = 'off';
ylim([0 380]);
xlim([-3 9]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[-3:9],'XTickLabel',{'','','','Arrival','','','Symptom onset','','','','First transmission event','',''},'Fontsize',18);
plot(linspace(0,MLET,2),[360 360], 'k','LineWidth',2); 
text(mean([0 MLET]), 380, [num2str(round(MLET,1)) ,' days'],'Fontsize',16,'HorizontalAlignment','center');
plot([0 0],[355 365], 'k','LineWidth',2); 
plot([MLET MLET],[355 365], 'k','LineWidth',2); 
plot(linspace(MLET,MLET+MLE,2),[360 360], 'k','LineWidth',2); 
plot([MLET+MLE MLET+MLE],[355 365], 'k','LineWidth',2); 
plot([MLET MLET],[355 365], 'k','LineWidth',2); 
text(mean([MLET MLET+MLE]), 380, [num2str(round(MLE,1)) ,' days'],'Fontsize',16,'HorizontalAlignment','center');
load('ArrivalToSymptomOnset.mat');
C=unique(ATtoSO);
mm=0;
for ii=1:length(C)   
    f=find(ATtoSO==C(ii)); 
    if(length(f)>mm)
        mm=length(f);
    end
end
 dy=linspace(87.500,262.5000,mm);
for ii=1:length(C)
    f=find(ATtoSO==C(ii)); 
    if(C(ii)<=0)
        scatter(C(ii).*ones(length(f),1),dy(1:length(f)),40,CC(1,:),'filled');
    else
        scatter(C(ii).*ones(length(f),1),dy(1:length(f)),40,CIP(1,:),'filled');
    end
end
text(-4.53,max(ylim)*1.1,'D','Fontsize',32,'FontWeight','bold');
text(-2,100,'A','Fontsize',32,'FontWeight','bold');