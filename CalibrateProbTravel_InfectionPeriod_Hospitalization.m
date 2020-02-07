%% Used to calibrate the probability of travel for the simulation based on hosptilaiztion
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

TT(:,(length(T)+1):(length(T)+length(TW)))=min(TT(:,(length(T)+1):(length(T)+length(TW))),INDX); % Restrict time based on the min of time of first medical and travel ban in wuhan (Do this because 
TT(:,(length(T)+length(TW)+1):end)=min(TT(:,(length(T)+length(TW)+1):end),INDX2); % Restrict time based on the min of time of first medical and travel ban in Hubei
%Used for the tima after symptom onset
TAO=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
TAOT=TAO; %Used in indexing the time of onset (could have just used TT to save memory...)
TAO(TAOT>=INDXMV)=TT(TAOT>=INDXMV)+TimeHospitalJan1(length(TAOT(TAOT>=INDXMV)));  % for symptom onset on jan1 or past sample from that distribtuion
TAO(TAOT<INDXMV)=TT(TAOT<INDXMV)+TimeHospitalDec31(length(TAOT(TAOT<INDXMV))); % for symptom onset before jan1 or past sample from that distribtuion
TBNS=TAO; % variable for integrating time ban 
TBNS(:,(length(T)+1):(length(T)+length(TW)))=min(TBNS(:,(length(T)+1):(length(T)+length(TW))),INDX); % Restrict time based on the min of time of first medical and travel ban in wuhan
TBNS(:,(length(T)+length(TW)+1):end)=min(TBNS(:,(length(T)+length(TW)+1):end),INDX2); % Restrict time based on the min of time of first medical and travel ban in Hubei

pc=0.01.*[0.3:0.01:3]; % Range used to calibrate the probability of travel
F=zeros(length(pc),1);
parfor mm=1:length(pc)
    ptravel=pc(mm); % set the travel probability
    UxT=zeros(NS1*NS2,maxE);

    D=(TT-E); % Calcualte the duration of the incubation period
    D(D<0)=0; % This ensures exposed people cannot leave during travel ban

    PItemp=wtc.*(1-(1-ptravel).^D); % probability of traveling outside of China
    
    TempT=[T TW TF]; % Do not need to truncate here as we need this for the time of symptom onset the probabailty takes care of the travel restriction
    for ii=1:maxE
       f=find(TempT==ii);
       if(~isempty(f))
            UxT(:,ii)=sum(PItemp(:,f),2); % sum those who have traveled and are now symtpmatic
       end
       for zz=1:NS2
          f=find(TBNS(zz,:)>ii); % find those who are symptomatic
          g=find(TempT(f)<=ii); % of those who still are are they in their symptomatic period
          dt=wtc.*ptravel*(1-ptravel).^(ii-TempT(f(g))); % probability
          UxT(zz,ii)=UxT(zz,ii)+sum(dt); % expected value
       end
    end    
    M=poisspdf(IncOutside',mean(UxT,1)); % calculate the likelihood at each point
    F(mm)=sum(log(M)); %sum the log likelihood
end
plot(pc,F)
save('Probability_Travel_Infection_Hospt.mat','F','pc');
