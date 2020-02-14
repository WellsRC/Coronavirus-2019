close all;
clear;
CC=[hex2rgb('#1E1F26');hex2rgb('#4D648D');hex2rgb('#D0E1F9');hex2rgb('#CB0000')];
startDateofSim = datenum('12-06-2019');% Start date
XTL=datestr([startDateofSim+[0:(68-1)]],'mm-dd-yy');
[IncC,IncW,IncH,IncO]=IncidenceData;
figure('units','normalized','outerposition',[0 0 1 1]);

subplot('Position',[0.065651260504202,0.64,0.921743697478991,0.3]); 
b=bar([1:68],[ IncW(:,2) IncH(:,2) IncC(:,2) IncO(:,2)],'stacked','LineStyle','none'); hold on;
for ii=1:4
        b(ii).FaceColor = 'flat';
        b(ii).CData = CC(ii,:);
end
box off;
% Index times of Key events
%Lockdown Wuhan
INDX=datenum('01-23-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero 
% Lockdown Hubei
INDX2=datenum('01-25-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
set(gca,'LineWidth',2,'tickdir','out','Fontsize',14,'XTick',[1:68],'XTickLabel',XTL,'Xminortick','off','Yminortick','on');
xtickangle(45);
xlim([1 68+0.5])
ylim([0 3500]);
yh=ylabel({'Incidence'},'Fontsize',18);
xlabel({'Date of symptom onset'},'Fontsize',18);
%Plot travel bans
plot([INDX INDX],[0 3500],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
plot([INDX2 INDX2],[0 3500],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
legend(b,{'Wuhan','Hubei (Outside Wuhan)','China (Outside Wuhan and Hubei)','International (Known Travel to China)'},'Location','NorthWest')
legend boxoff;
text(yh.Extent(1),max(ylim)*1.1,'A','Fontsize',32,'FontWeight','bold');
subplot('Position',[0.065651260504202,0.165,0.921743697478991,0.3]); 
load('ReportedIncidenceChina-Dec62019=t0');
load('ReportedIncidenceOther-WHO-Dec62019=t0');
load('ReportedIncidenceWuhan-Dec62019=t0');
load('ReportedIncidenceHubei-Dec62019=t0');

b=bar([1:68],[ RWuhan(:,2) RHubei(:,2) RChina(:,2) OWHORep(:,2)],'stacked','LineStyle','none'); hold on;
for ii=1:4
        b(ii).FaceColor = 'flat';
        b(ii).CData = CC(ii,:);
end
box off;
% Index times of Key events
%Lockdown Wuhan
INDX=datenum('01-23-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero 
% Lockdown Hubei
INDX2=datenum('01-25-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
set(gca,'LineWidth',2,'tickdir','out','Fontsize',14,'XTick',[1:68],'XTickLabel',XTL,'Xminortick','off','Yminortick','on','YTick',[0:500:4500]);
xtickangle(45);
xlim([1 68+0.5])
ylim([0 4500]);
ylabel({'Incidence'},'Fontsize',18);
xlabel({'Date of report'},'Fontsize',18);
%Plot travel bans
plot([INDX INDX],[0 4500],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
plot([INDX2 INDX2],[0 4500],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
legend(b,{'Wuhan','Hubei (Outside Wuhan)','China (Outside Wuhan and Hubei)','International (Known Travel to China)'},'Location','NorthWest')
legend boxoff;

text(yh.Extent(1),max(ylim)*1.1,'B','Fontsize',32,'FontWeight','bold');