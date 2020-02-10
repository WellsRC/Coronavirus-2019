%% Geenrates Figure S1
clear;
% Determine the weight for flgihts outside of China
load('Weight_Flights.mat','FlightAll')
tf=strcmp({'China'},{FlightAll{:,1}});
wtc=1-[FlightAll{tf,2}]; % the weight for flights outside of China
% Load the incidence data
[IncC,IncW,IncH,IncO]=IncidenceData;
% Number of samples to generate
NS1=1; % NS1 set to one as concerned with mle of incubatino period
NS2=10^3; 

% Time vectors for symptoms onset
T=[]; % Outside wuhan and Hubei
TW=[]; % In Wuhan
TF=[]; % In Hubei
% Set the time of symptom onset for each case (Note: Day index 1 correspnds to Dec 6)
%China (outside Wuhan and Hubei)
for ii=1:length(IncC(:,1))
   T=[T IncC(ii,1).*ones(1,IncC(ii,2))]; 
end
% International cases
for ii=1:length(IncO(:,1))
   T=[T IncO(ii,1).*ones(1,IncO(ii,2))]; 
end
%Wuhan
for ii=1:length(IncW(:,1))
   TW=[TW IncW(ii,1).*ones(1,IncW(ii,2))]; 
end
%Hubeia
for ii=1:length(IncH(:,1))
   TF=[TF IncH(ii,1).*ones(1,IncH(ii,2))]; 
end
% The daily incidence for international cases
IncOutside=IncO(:,2);

% Index times of Key events
%Lockdown Wuhan
INDX=datenum('01-23-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero 
% Lockdown Hubei
INDX2=datenum('01-25-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
%Switch in in distribution for time to seek medical
INDXMV=datenum('01-1-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero

% Minimum day index (as we sample back for Exposed period
minE=-22;%min(E(:));
% Maximim index
maxE=max(IncO(:,1));

% Sample duration of the incubation period
IP=zeros(NS1*NS2,length(T)+length(TW)+length(TF));
for nn=1:NS1    
    for ii=1:NS2
        [IP(ii+NS2.*(nn-1),:),~] = IncubationDist(5.2,length(T)+length(TW)+length(TF)); % Samples from rounded log-normal with a mean of 5.2
    end
end
% Time of infection
E=repmat([T TW TF],NS2,1)-IP; 
% Time of symptom onset
%(this preamble was copied were I run all cases serial in the simulatino later)
TT=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];

TT(:,(length(T)+1):(length(T)+length(TW)))=min(TT(:,(length(T)+1):(length(T)+length(TW))),INDX); % Restrict time based on the min of time of first medical and travel ban in wuhan
TT(:,(length(T)+length(TW)+1):end)=min(TT(:,(length(T)+length(TW)+1):end),INDX2); % Restrict time based on the min of time of first medical and travel ban in Hubei
%Used for the tima after symptom onset
TAO=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
TAOT=TAO; %Used in indexing the time of onset (could have just used TT to save memory...)
TAO(TAOT>=INDXMV)=TT(TAOT>=INDXMV)+TimeMedJan1(length(TAOT(TAOT>=INDXMV)));  % for symptom onset on jan1 or past sample from that distribtuion
TAO(TAOT<INDXMV)=TT(TAOT<INDXMV)+TimeMedDec31(length(TAOT(TAOT<INDXMV))); % for symptom onset before jan1 or past sample from that distribtuion
TBNS=TAO; % variable for integrating time ban 
TBNS(:,(length(T)+1):(length(T)+length(TW)))=min(TBNS(:,(length(T)+1):(length(T)+length(TW))),INDX); % Restrict time based on the min of time of first medical and travel ban in wuhan
TBNS(:,(length(T)+length(TW)+1):end)=min(TBNS(:,(length(T)+length(TW)+1):end),INDX2); % Restrict time based on the min of time of first medical and travel ban in Hubei

% Load calibrate dprobability
load('Probability_Travel_Infection.mat','F','pc');
w=exp(F)./sum(exp(F));
wc=cumsum(w);

pmb=[pc(F==max(F)) 0 0];
r=rand(10^4,1);
spc=zeros(10^4,1);

for ii=1:10^4
    f=find(r(ii)<=wc);
    f=f(1);
    spc(ii)=pc(f);
end

pmb(2)=prctile(spc,2.5);
pmb(3)=prctile(spc,97.5);


F=zeros(3,length([1:maxE])); % only need to go from one to maxE as no infectious case is past one and we are looking at the symptomatic cases only
for mm=1:length(pmb)
    ptravel=pmb(mm);
    UxT=zeros(NS1*NS2,maxE);

    D=(TT-E);
    D(D<0)=0;

    PItemp=wtc.*(1-(1-ptravel).^D); % Travel before symptoms appear
    
    TempT=[T TW TF];
    for ii=1:maxE
       f=find(TempT==ii);
       if(~isempty(f))
            UxT(:,ii)=sum(PItemp(:,f),2); % time of sympaots based on travel over incubation period
       end
       for zz=1:NS2 % travel during symptoms
          f=find(TBNS(zz,:)>ii);
          g=find(TempT(f)<=ii);
          dt=wtc.*ptravel*(1-ptravel).^(ii-TempT(f(g)));
          UxT(zz,ii)=UxT(zz,ii)+sum(dt);
       end
    end    
    F(mm,:)=mean(UxT,1);
end
startDateofSim = datenum('12-06-2019');% Start date
XTL=datestr([startDateofSim+[0:4:(maxE-1)]],'mm-dd-yy');

figure('units','normalized','outerposition',[0 0 1 1]);
% Plot the patches
for ii=1:maxE
   patch(ii+[-0.35 0.35 0.35 -0.35], [F(2,ii) F(2,ii) F(3,ii) F(3,ii)],'k','LineStyle','none','Facealpha',0.3);hold on
   plot(ii+linspace(-0.35,0.35,2),F(1,ii).*ones(1,2),'k','LineWidth',2);    
   scatter(IncO(ii,1),IncO(ii,2),40,'r','filled');
end
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:4:maxE],'XTickLabel',XTL,'Xminortick','on','Yminortick','on');
xtickangle(45);
xlim([1 maxE+0.5])
ylim([0 30]);
ylabel({'Incidence'},'Fontsize',18);
xlabel({'Date'},'Fontsize',18);
%Plot travel bans
plot([INDX INDX],[0 30],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
plot([INDX2 INDX2],[0 30],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);