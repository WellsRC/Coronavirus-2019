function [IncCF,IncWF,IncHF,IncOF]=IncidenceData
ENDINDX=datenum('01-31-2020')-datenum('12-6-2019'); % If we had a data point for Dec 6 it would have the time stamp of day zero, where Dec 8 would have a week index of two
IncCF=zeros(ENDINDX+1,2); % China 
IncWF=zeros(ENDINDX+1,2); % Wuhan
IncHF=zeros(ENDINDX+1,2); % Hubei
IncOF=zeros(ENDINDX+1,2); % Other
IncCF(:,1)=[1:ENDINDX+1];
IncWF(:,1)=[1:ENDINDX+1];
IncHF(:,1)=[1:ENDINDX+1];
IncOF(:,1)=[1:ENDINDX+1];
%% Load time of symptom onset
load('DailyIncidenceWuhan-NEJM-Dec62019=t0');
IncWF(WNEJM(:,1),2)=WNEJM(:,2); % Need to add one as the week index for Dec 6 would be zero
load('DailyIncidenceOther-WHO-Dec62019=t0');
IncOF(OWHO(:,1),2)=OWHO(:,2);
%% Reported time of incidence (Adjust backwards based on time of first medical visit since later in the outbreak)
% Seed random number generator for the results remain consistent
rng(20200130);
load('DailyReportedInidenceChina-Affan-Dec62019=t0');
for ii=1:length(RChinaAffan(:,1))
    [TH]=TimeMedJan1(RChinaAffan(ii,2));
    for jj=1:length(TH)
        IncCF(RChinaAffan(ii,1)-TH(jj),2)=IncCF(RChinaAffan(ii,1)-TH(jj),2)+1;
    end
end

load('DailyReportedInidenceOther-Affan-Dec62019=t0');
for ii=1:length(ROtherAffan(:,1))
    [TH]=TimeMedJan1(ROtherAffan(ii,2));
    for jj=1:length(TH)
        IncOF(ROtherAffan(ii,1)-TH(jj),2)=IncOF(ROtherAffan(ii,1)-TH(jj),2)+1;
    end
end

load('DailyReportedInidenceWuhan-Affan-Dec62019=t0');
for ii=1:length(RWAffan(:,1))
    [TH]=TimeMedJan1(RWAffan(ii,2));
    for jj=1:length(TH)
        IncWF(RWAffan(ii,1)-TH(jj),2)=IncWF(RWAffan(ii,1)-TH(jj),2)+1;
    end
end


load('DailyReportedInidenceHubei-Affan-Dec62019=t0');
for ii=1:length(RHAffan(:,1))
    [TH]=TimeMedJan1(RHAffan(ii,2));
    for jj=1:length(TH)
        IncHF(RHAffan(ii,1)-TH(jj),2)=IncHF(RHAffan(ii,1)-TH(jj),2)+1;
    end
end
end