% Runs analysis for Figure 1
clear;

% pobj=parpool(20);
% Determine the weight for flgihts outside of China
load('Weight_Flights.mat','FlightAll','Flight_NW')
tf=strcmp({'China'},{FlightAll{:,1}});
wtc1=1-[FlightAll{tf,2}]; % the weight for flights outside of China
tf=strcmp({'China'},{Flight_NW{:,1}});
wtc2=1-[Flight_NW{tf,2}]; % the weight for flights outside of China
% Load the incidence data
[IncC,IncW,IncH,IncO]=IncidenceData;
% Number of samples to generate
NS1=10^3; 
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

% Estimates paramter for gamma with mean 5.2 
options=optimset('MaxIter',10^6,'TolFun',10^(-16),'TolX',10^(-16),'display','off');
[gp]=fmincon(@(x)((gaminv(0.025,x(1),5.2./x(1))-4.1).^2+(gaminv(0.975,x(1),5.2./x(1))-7).^2),100,[],[],[],[],0,1000,[],options);
mun=gamrnd(gp,5.2/gp,NS1,1);

% loads calibrated prob of travel
load('Probability_Travel_Hospital_6733.mat','F','pc');
w=exp(F)./sum(exp(F)); % Calcualted likelihood wieght
wc=cumsum(w); % cumualtive sum for smapling

ptravel=pc(F==max(F)); % mle


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


%% MLE Estimates
% Initialize for reuslts 
UxT=zeros(NS2,length([minE:maxE]));
UxS=zeros(NS2,length([minE:maxE]));
UxNS=zeros(NS2,length([minE:maxE]));
UxTNS=zeros(NS2,length([minE:maxE]));

PxT=ones(NS2,length([minE:maxE]));
PxS=ones(NS2,length([minE:maxE]));
CPxTS=ones(NS2,length([minE:maxE]));
CPxTNS=ones(NS2,length([minE:maxE]));
PxNS=ones(NS2,length([minE:maxE]));
PxTNS=ones(NS2,length([minE:maxE]));

% Sample duration of the incubation period
IP=zeros(NS2,length(T)+length(TW)+length(TF));
for ii=1:NS2
    [IP(ii,:),~] = IncubationDist(5.2,length(T)+length(TW)+length(TF)); % Samples from rounded log-normal with a mean of 5.2
end


% Travel Ban and screening
E=repmat([T TW TF],NS2,1)-IP; 
TT=[repmat([T],NS2,1) min(repmat([TW],NS2,1),INDX) min(repmat([TF],NS2,1),INDX2)]; % Not subtracting one as when calcuating below we look at days exclusive before and do not include this day
TNR=repmat([T TW TF],NS2,1);
TAO=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
TAOT=TAO;
TAO(TAOT>=INDXMV)=TNR(TAOT>=INDXMV)+TimeHospitalJan1(length(TAOT(TAOT>=INDXMV))); 
TAO(TAOT<INDXMV)=TNR(TAOT<INDXMV)+TimeHospitalDec31(length(TAOT(TAOT<INDXMV))); 
TBNS=TAO;
TBNS(:,(length(T)+1):(length(T)+length(TW)))=min(TBNS(:,(length(T)+1):(length(T)+length(TW))),INDX);
TBNS(:,(length(T)+length(TW)+1):end)=min(TBNS(:,(length(T)+length(TW)+1):end),INDX2);
TempZ=zeros(size(E));
TempZ2=zeros(size(E));
for ii=minE:maxE
    if(ii>=INDX)
      wtc=wtc2; 
    else          
      wtc=wtc1; 
    end
   for mm=1:NS2          
       f=find(E(mm,:)<=ii);
      g=find(TT(mm,f)>ii);
      dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
      UxT(mm,ii-(minE)+1)=sum(dt);
      PxT(mm,ii-(minE)+1)=(1-prod(1-dt));
      TempZ(mm,f(g))=TempZ(mm,f(g))+dt;
      CPxTS(mm,ii-(minE)+1)=(1-prod(1-TempZ(mm,:)));


      g=find(TNR(mm,f)>ii);
      dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
      UxS(mm,ii-(minE)+1)=sum(dt);
      PxS(mm,ii-(minE)+1)=(1-prod(1-dt));

      g=find(TAO(mm,f)>ii);
      dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
      UxNS(mm,ii-(minE)+1)=sum(dt);
      PxNS(mm,ii-(minE)+1)=(1-prod(1-dt));

      g=find(TBNS(mm,f)>ii);
      dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
      UxTNS(mm,ii-(minE)+1)=sum(dt);
      PxTNS(mm,ii-(minE)+1)=(1-prod(1-dt));
      TempZ2(mm,f(g))=TempZ2(mm,f(g))+dt;
      CPxTNS(mm,ii-(minE)+1)=(1-prod(1-TempZ2(mm,:)));

   end
end    

MLExTS=mean(UxT,1);
MLExS=mean(UxS,1);
MLExNS=mean(UxNS,1);
MLExTNS=mean(UxTNS,1);

MPTS=mean(PxT,1);
MPS=mean(PxS,1);
MPNS=mean(PxNS,1);
MPTNS=mean(PxTNS,1);

MCPTS=mean(CPxTS,1);
MCPTNS=mean(CPxTNS,1);
save('Daily_Prob_Expect_Hospital_6733.mat','MLExTS','MLExS','MLExNS','MLExTNS','MPTS','MPS','MPNS','MPTNS','MCPTS','MCPTNS'););
%% Uncertainty Estimates

UMLExTS=zeros(NS1,length(MLExTS));
UMLExS=zeros(NS1,length(MLExTS));
UMLExNS=zeros(NS1,length(MLExTS));
UMLExTNS=zeros(NS1,length(MLExTS));


UMCPTS=zeros(NS1,length(MLExTS));
UMCPTNS=zeros(NS1,length(MLExTS));

UMPTS=zeros(NS1,length(MLExTS));
UMPS=zeros(NS1,length(MLExTS));
UMPNS=zeros(NS1,length(MLExTS));
UMPTNS=zeros(NS1,length(MLExTS));
r=rand(NS1,1);
spc=zeros(NS1,1);

for ii=1:NS1
    f=find(r(ii)<=wc);
    f=f(1);
    spc(ii)=pc(f);
end

parfor ss=1:NS1
    ptravel=spc(ss);
    UxT=zeros(NS2,length([minE:maxE]));
    UxS=zeros(NS2,length([minE:maxE]));
    UxNS=zeros(NS2,length([minE:maxE]));
    UxTNS=zeros(NS2,length([minE:maxE]));
    IP=zeros(NS2,length(T)+length(TW)+length(TF));

    PxT=ones(NS2,length([minE:maxE]));
    PxS=ones(NS2,length([minE:maxE]));
    CPxTS=ones(NS2,length([minE:maxE]));
    CPxTNS=ones(NS2,length([minE:maxE]));
    PxNS=ones(NS2,length([minE:maxE]));
    PxTNS=ones(NS2,length([minE:maxE]));

    % Sample incubation period
    for ii=1:NS2
        [IP(ii,:),~] = IncubationDist(mun(ss),length(T)+length(TW)+length(TF));
    end

    % Travel Ban and screening
    E=repmat([T TW TF],NS2,1)-IP; 
   TT=[repmat([T],NS2,1) min(repmat([TW],NS2,1),INDX) min(repmat([TF],NS2,1),INDX2)]; % Not subtracting one as when calcuating below we look at days exclusive before and do not include this day
TNR=repmat([T TW TF],NS2,1);
TAO=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
TAOT=TAO;
TAO(TAOT>=INDXMV)=TNR(TAOT>=INDXMV)+TimeHospitalJan1(length(TAOT(TAOT>=INDXMV))); 
TAO(TAOT<INDXMV)=TNR(TAOT<INDXMV)+TimeHospitalDec31(length(TAOT(TAOT<INDXMV))); 
TBNS=TAO;
TBNS(:,(length(T)+1):(length(T)+length(TW)))=min(TBNS(:,(length(T)+1):(length(T)+length(TW))),INDX);
TBNS(:,(length(T)+length(TW)+1):end)=min(TBNS(:,(length(T)+length(TW)+1):end),INDX2);
TempZ=zeros(size(E));
TempZ2=zeros(size(E));
    for ii=minE:maxE
        if(ii>=INDX)
          wtc=wtc2; 
        else          
          wtc=wtc1; 
        end
       for mm=1:NS2         
           f=find(E(mm,:)<=ii);
          g=find(TT(mm,f)>ii);
          dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
          UxT(mm,ii-(minE)+1)=sum(dt);
          PxT(mm,ii-(minE)+1)=(1-prod(1-dt));
          TempZ(mm,f(g))=TempZ(mm,f(g))+dt;
          CPxTS(mm,ii-(minE)+1)=(1-prod(1-TempZ(mm,:)));


          g=find(TNR(mm,f)>ii);
          dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
          UxS(mm,ii-(minE)+1)=sum(dt);
          PxS(mm,ii-(minE)+1)=(1-prod(1-dt));

          g=find(TAO(mm,f)>ii);
          dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
          UxNS(mm,ii-(minE)+1)=sum(dt);
          PxNS(mm,ii-(minE)+1)=(1-prod(1-dt));

          g=find(TBNS(mm,f)>ii);
          dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
          UxTNS(mm,ii-(minE)+1)=sum(dt);
          PxTNS(mm,ii-(minE)+1)=(1-prod(1-dt));
          TempZ2(mm,f(g))=TempZ2(mm,f(g))+dt;
          CPxTNS(mm,ii-(minE)+1)=(1-prod(1-TempZ2(mm,:)));

       end
    end  

    UMLExTS(ss,:)=mean(UxT,1);
    UMLExS(ss,:)=mean(UxS,1);
    UMLExNS(ss,:)=mean(UxNS,1);
    UMLExTNS(ss,:)=mean(UxTNS,1);
    UMPTS(ss,:)=mean(PxT,1);
    UMPS(ss,:)=mean(PxS,1);
    UMPNS(ss,:)=mean(PxNS,1);
    UMPTNS(ss,:)=mean(PxTNS,1);
    UMCPTS(ss,:)=mean(CPxTS,1);
    UMCPTNS(ss,:)=mean(CPxTNS,1);
end
save('Daily_Prob_Expect_Hospital_6733.mat','UMLExTS','UMLExS','UMLExNS','UMLExTNS','UMPTS','UMPS','UMPNS','UMPTNS','UMCPTS','UMCPTNS','MLExTS','MLExS','MLExNS','MLExTNS','MPTS','MPS','MPNS','MPTNS','MCPTS','MCPTNS'););
clear;