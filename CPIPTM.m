%% Used to calibrate the probability of travel for the simulation (Testing the effect of inclusing the reporting date)
clear;

pobj=parpool(10);
% Determine the weight for flgihts outside of China
load('Weight_Flights.mat')
tf=strcmp({'China'},{FlightAll{:,1}});
w1=1-[FlightAll{tf,2}]; % the weight for flights outside of China

tf=strcmp({'China'},{Flight_NW{:,1}});
w2=1-[Flight_NW{tf,2}]; % the weight for flights outside of China

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
load('ReportedTimeIncidence.mat');
RepI=RepI(:,2);
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
TAO(TAOT>=INDXMV)=TT(TAOT>=INDXMV)+TimeMedJan1(length(TAOT(TAOT>=INDXMV)));  % for symptom onset on jan1 or past sample from that distribtuion (Note: here it is not an issue we truncated TT indexes before as we truncate this later as well and we are adding on to of the time of the travle ban)
TAO(TAOT<INDXMV)=TT(TAOT<INDXMV)+TimeMedDec31(length(TAOT(TAOT<INDXMV))); % for symptom onset before jan1 or past sample from that distribtuion (Note: here it is not an issue we truncated TT indexes before as we truncate this later as well and we are adding on to of the time of the travle ban)
TBNS=TAO; % variable for integrating time ban 
TBNS(:,(length(T)+1):(length(T)+length(TW)))=min(TBNS(:,(length(T)+1):(length(T)+length(TW))),INDX); % Restrict time based on the min of time of first medical and travel ban in wuhan
TBNS(:,(length(T)+length(TW)+1):end)=min(TBNS(:,(length(T)+length(TW)+1):end),INDX2); % Restrict time based on the min of time of first medical and travel ban in Hubei

pc=0.01.*[0.3:0.01:3]; % Range used to calibrate the probability of travel
F=zeros(72-36+1,length(pc));
parfor mm=1:length(pc)
    ptravel=pc(mm); % set the travel probability
    UxT2=zeros(NS1*NS2,maxE);
    UxT=zeros(NS1*NS2,maxE);
    PItemp=zeros(size(E)); % probability of traveling outside of China
    for ii=minE:maxE
       if(ii>=INDX)
          wtc=w2; 
       else          
          wtc=w1; 
       end
       for jj=1:NS2
           f=find(TT(jj,:)>ii); % find those liekly in infected period (THIS TAKES CARE OF TRAVEL BAN)
          g=find(E(jj,f)<=ii); % find thos that are in infected period (WONT TAKE PEOPLE IN TRAVEL BAN e.g. if travel ban is ii the ywill not be selected in line above)
          dt=wtc.*ptravel*(1-ptravel).^(ii-E(jj,f(g))); % probability
          PItemp(jj,f(g))=PItemp(jj,f(g))+dt;
       end
    end
    
    TempT=[T TW TF]; % Do not need to truncate here as we need this for the time of symptom onset the probabailty takes care of the travel restriction a exposed people can travel before the travel ban and then show symptoms
    for ii=1:maxE
       if(ii>=INDX)
          wtc=w2; 
       else          
          wtc=w1; 
       end
       f=find(TempT==ii); % Find the times of symptom onset
       if(~isempty(f))
            UxT2(:,ii)=sum(PItemp(:,f),2); % sum those who have traveled and are now symtpmatic
            UxT(:,ii)=sum(PItemp(:,f),2); % sum those who have traveled and are now symtpmatic
       end
       for zz=1:NS2
          f=find(TBNS(zz,:)>ii); % find those who are symptomatic
          g=find(TempT(f)<=ii); % of those who still are are they in their symptomatic period
          dt=wtc.*ptravel*(1-ptravel).^(ii-TempT(f(g))); % probability
          UxT2(zz,ii)=UxT2(zz,ii)+sum(dt); % expected value
       end
    end    
    TM1=mean(UxT,1);
    TM2=mean(UxT2,1);
    for ttt=1:37
        M1=poisspdf(IncOutside(1:(35+ttt))',TM1(1:(35+ttt))); % calculate the likelihood at each point for date of onset for those who travel during incubation period
        M2=poisspdf(RepI(1:(35+ttt))',TM2(1:(35+ttt))); % calculate the likelihood at each point for reported data for those who travell during sympmatic period
        F(ttt,mm)=(20/30).*sum(log(M1))+(10/30).*sum(log(M2)); %sum the log likelihood (Weighted to account for those who travel during incubation period and those who travell during symptkmatic period)
    end
end
%plot(pc,F)
save('Probability_Travel_Infection_6733_Mulit.mat','F','pc');
