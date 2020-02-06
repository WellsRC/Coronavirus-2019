clear;
load('Weight_Flights.mat','FlightAll')
tf=strcmp({'China'},{FlightAll{:,1}});
wtc=1-[FlightAll{tf,2}]; % Not to china
%pobj=parpool(20);
[IncC,IncW,IncH,IncO]=IncidenceData;
NS1=10^3;
NS2=10^3;
T=[];
TW=[];
TF=[];
for ii=1:length(IncC(:,1))
   T=[T IncC(ii,1).*ones(1,IncC(ii,2))]; 
end
for ii=1:length(IncO(:,1))
   T=[T IncO(ii,1).*ones(1,IncO(ii,2))]; 
end

for ii=1:length(IncW(:,1))
   TW=[TW IncW(ii,1).*ones(1,IncW(ii,2))]; 
end


for ii=1:length(IncH(:,1))
   TF=[TF IncH(ii,1).*ones(1,IncH(ii,2))]; 
end
IncOutside=IncO(:,2);
options=optimset('MaxIter',10^6,'TolFun',10^(-16),'TolX',10^(-16),'display','off');
[gp]=fmincon(@(x)((gaminv(0.025,x(1),5.2./x(1))-4.1).^2+(gaminv(0.975,x(1),5.2./x(1))-7).^2),100,[],[],[],[],0,1000,[],options);
gaminv(0.025,gp,5.2./gp)
gaminv(0.975,gp,5.2./gp)
mun=gamrnd(gp,5.2/gp,NS1,1);


load('Probability_Travel_Infection.mat','F','pc');
w=exp(F)./sum(exp(F));
wc=cumsum(w);

ptravel=pc(F==max(F));


INDX=datenum('01-23-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
INDX2=datenum('01-25-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
INDXMV=datenum('01-1-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero

minE=-22;%min(E(:));
maxE=max(IncC(:,1));

TA=ones(1,length([minE:maxE]));
TAN=ones(1,length([minE:maxE]));
TA(INDX-minE:end)=0;

%% MLE Estimates

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
    [IP(ii,:),~] = IncubationDist(5.2,length(T)+length(TW)+length(TF));
end

% Travel Ban and screening
E=repmat([T TW TF],NS2,1)-IP; 
TT=[repmat([T],NS2,1) min(repmat([TW],NS2,1),INDX) min(repmat([TF],NS2,1),INDX2)]; % Not subtracting one as when calcuating below we look at days exclusive before and do not include this day
TNR=repmat([T TW TF],NS2,1);
TAO=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
TAOT=TAO;
TAO(TAOT>=INDXMV)=TNR(TAOT>=INDXMV)+TimeMedJan1(length(TAOT(TAOT>=INDXMV))); 
TAO(TAOT<INDXMV)=TNR(TAOT<INDXMV)+TimeMedDec31(length(TAOT(TAOT<INDXMV))); 
TBNS=TAO;
TBNS(:,(length(T)+1):(length(T)+length(TW)))=min(TBNS(:,(length(T)+1):(length(T)+length(TW))),INDX);
TBNS(:,(length(T)+length(TW)+1):end)=min(TBNS(:,(length(T)+length(TW)+1):end),INDX2);
for ii=minE:maxE
   for mm=1:NS2
      f=find(TT(mm,:)>ii);
      g=find(E(mm,f)<=ii);
      dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
      UxT(mm,ii-(minE)+1)=UxT(mm,ii-(minE)+1)+sum(dt);
      PxT(mm,ii-(minE)+1)=PxT(mm,ii-(minE)+1).*(1-prod(1-dt));
      gg=find(E(mm,:)<=ii);
      temps=min(TT(mm,gg),ii+1);
      dts=(1-prod(1-wtc.*(1-(1-ptravel).^temps)));
      CPxTS(mm,ii-(minE)+1)=CPxTS(mm,ii-(minE)+1).*dts;
      
      
      f=find(TNR(mm,:)>ii);
      g=find(E(mm,f)<=ii);
      dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
      UxS(mm,ii-(minE)+1)=UxS(mm,ii-(minE)+1)+sum(dt);
      PxS(mm,ii-(minE)+1)=PxS(mm,ii-(minE)+1).*(1-prod(1-dt));
      
      f=find(TAO(mm,:)>ii);
      g=find(E(mm,f)<=ii);
      dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
      UxNS(mm,ii-(minE)+1)=UxNS(mm,ii-(minE)+1)+sum(dt);
      PxNS(mm,ii-(minE)+1)=PxNS(mm,ii-(minE)+1).*(1-prod(1-dt));
      
      f=find(TBNS(mm,:)>ii);
      g=find(E(mm,f)<=ii);
      dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
      UxTNS(mm,ii-(minE)+1)=UxTNS(mm,ii-(minE)+1)+sum(dt);
      PxTNS(mm,ii-(minE)+1)=PxTNS(mm,ii-(minE)+1).*(1-prod(1-dt));
      gg=find(E(mm,:)<=ii);
      temps=min(TBNS(mm,gg),ii+1);
      dts=(1-prod(1-wtc.*(1-(1-ptravel).^temps)));
      CPxTNS(mm,ii-(minE)+1)=CPxTNS(mm,ii-(minE)+1).*dts;
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
save('Daily_Prob_Expect.mat');
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
    TAO(TAOT>=INDXMV)=TNR(TAOT>=INDXMV)+TimeMedJan1(length(TAOT(TAOT>=INDXMV))); 
    TAO(TAOT<INDXMV)=TNR(TAOT<INDXMV)+TimeMedDec31(length(TAOT(TAOT<INDXMV))); 
    TBNS=TAO;
    TBNS(:,(length(T)+1):(length(T)+length(TW)))=min(TBNS(:,(length(T)+1):(length(T)+length(TW))),INDX);
    TBNS(:,(length(T)+length(TW)+1):end)=min(TBNS(:,(length(T)+length(TW)+1):end),INDX2);
    for ii=minE:maxE
       for mm=1:NS2
          f=find(TT(mm,:)>ii);
          g=find(E(mm,f)<=ii);
          dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
          UxT(mm,ii-(minE)+1)=UxT(mm,ii-(minE)+1)+sum(dt);
          PxT(mm,ii-(minE)+1)=PxT(mm,ii-(minE)+1).*(1-prod(1-dt));
          gg=find(E(mm,:)<=ii);
          temps=min(TT(mm,gg),ii+1);
          dts=(1-prod(1-wtc.*(1-(1-ptravel).^temps)));
          CPxTS(mm,ii-(minE)+1)=CPxTS(mm,ii-(minE)+1).*dts;


          f=find(TNR(mm,:)>ii);
          g=find(E(mm,f)<=ii);
          dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
          UxS(mm,ii-(minE)+1)=UxS(mm,ii-(minE)+1)+sum(dt);
          PxS(mm,ii-(minE)+1)=PxS(mm,ii-(minE)+1).*(1-prod(1-dt));

          f=find(TAO(mm,:)>ii);
          g=find(E(mm,f)<=ii);
          dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
          UxNS(mm,ii-(minE)+1)=UxNS(mm,ii-(minE)+1)+sum(dt);
          PxNS(mm,ii-(minE)+1)=PxNS(mm,ii-(minE)+1).*(1-prod(1-dt));

          f=find(TBNS(mm,:)>ii);
          g=find(E(mm,f)<=ii);
          dt=wtc.*ptravel*(1-ptravel).^(ii-E(mm,f(g)));
          UxTNS(mm,ii-(minE)+1)=UxTNS(mm,ii-(minE)+1)+sum(dt);
          PxTNS(mm,ii-(minE)+1)=PxTNS(mm,ii-(minE)+1).*(1-prod(1-dt));
          gg=find(E(mm,:)<=ii);
          temps=min(TBNS(mm,gg),ii+1);
          dts=(1-prod(1-wtc.*(1-(1-ptravel).^temps)));
          CPxTNS(mm,ii-(minE)+1)=CPxTNS(mm,ii-(minE)+1).*dts;
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
save('Daily_Prob_Expect.mat');
clear;