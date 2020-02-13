%% Calcualtes the Probabilities for the specific countries
clear;
% Determine the weight for flgihts outside of China


% Load the incidence data
[IncC,IncW,IncH,IncO]=IncidenceData;
% Number of samples to generate
NS1=1; 
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

% Load the probability of travel that was calibrated
load('Probability_Travel_Infection.mat','F','pc');
ptravel=pc(F==max(F)); % set to the mle

% Index times of Key events
%Lockdown Wuhan
INDX=datenum('01-23-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero 
% Lockdown Hubei
INDX2=datenum('01-25-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
%Switch in in distribution for time to seek medical
INDXMV=datenum('01-1-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero


minE=-22;%min(E(:));
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
TT=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)]; % Do not need to truncate here as we go sted by step and use the travel ban
%Used for the tima after symptom onset
TAO=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
TAOT=TAO; %Used in indexing the time of onset (could have just used TT to save memory...)
TAO(TAOT>=INDXMV)=TT(TAOT>=INDXMV)+TimeMedJan1(length(TAOT(TAOT>=INDXMV)));  % for symptom onset on jan1 or past sample from that distribtuion
TAO(TAOT<INDXMV)=TT(TAOT<INDXMV)+TimeMedDec31(length(TAOT(TAOT<INDXMV))); % for symptom onset before jan1 or past sample from that distribtuion
TBNS=TAO; % variable for integrating time ban 
TBNS(:,(length(T)+1):(length(T)+length(TW)))=min(TBNS(:,(length(T)+1):(length(T)+length(TW))),INDX); % Restrict time based on the min of time of first medical and travel ban in wuhan
TBNS(:,(length(T)+length(TW)+1):end)=min(TBNS(:,(length(T)+length(TW)+1):end),INDX2); % Restrict time based on the min of time of first medical and travel ban in Hubei

% Determine the flight weights of the conuntries we need
load('Weight_Flights.mat')
w=zeros(length([FC{:,2}]),1);
for ii=1:length(w)    
    tf = strcmp({FC{ii,1}},{FlightAll{:,1}}); % find country specific name
    w(ii)=FlightAll{tf,2}; % define the wieght
end
w2=zeros(length([FC{:,2}]),1);
for ii=1:length(w2)    
    tf = strcmp({FC{ii,1}},{Flight_NW{:,1}}); % find country specific name
    w2(ii)=Flight_NW{tf,2}; % define the wieght
end

MLE=zeros(length(w),length([minE:maxE])); % Daily probability 
MLEP=zeros(length(w),length([minE:maxE])); % Cumualtive probability
parfor ww=1:length([FC{:,2}])
    
    TTT=zeros(NS2,length([T TW TF]));
    CPxTNS=ones(NS2,length([minE:maxE])); % initialize cumulative prob.
    PxTNS=ones(NS2,length([minE:maxE])); % initialize daily prob.
    for ii=minE:maxE
       if(ii>=INDX)
          wtc=w2(ww); 
       else          
          wtc=w(ww); 
       end
       for mm=1:NS2
           f=find(TBNS(mm,:)>ii); % find those liekly in infected period (THIS TAKES CARE OF TRAVEL BAN)
          g=find(E(mm,f)<=ii); % find thos that are in infected period (WONT TAKE PEOPLE IN TRAVEL BAN e.g. if travel ban is ii the ywill not be selected in line above)
          dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g))); % calcualte the probability of each case
          TTT(mm,f(g))=TTT(mm,f(g))+dt;
          PxTNS(mm,ii-(minE)+1)=(1-prod(1-dt)); % calaculte total probabaility 
          dts=(1-prod(1-TTT(mm,:))); % calcualte probability for eahc
          CPxTNS(mm,ii-(minE)+1)=dts; % put to cumualtive
       end
    end
    MLE(ww,:)=mean(CPxTNS,1);
    MLEP(ww,:)=mean(PxTNS,1);
end

    save(['Weighted_Travel_Infectious_Country.mat'],'MLE','w','w2','MLEP');