clear;
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

ptravelB=pc(F==max(F));
INDX=datenum('01-23-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
INDX2=datenum('01-25-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
INDXMV=datenum('01-1-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero


minE=-22;%min(E(:));
maxE=max(IncO(:,1));

%% MLE Estimates
IP=zeros(NS2,length(T)+length(TW)+length(TF));


% Sample incubation period
for ii=1:NS2
    [IP(ii,:),~] = IncubationDist(5.2,length(T)+length(TW)+length(TF));
end

% Travel Ban Infectious
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

load('Weight_Flights.mat')
w=zeros(length([FC{:,2}]),1);
for ii=1:length(w)    
    tf = strcmp({FC{ii,1}},{FlightAll{:,1}});
    w(ii)=FlightAll{tf,2};
end
w=unique(w);
MLE=zeros(length(w),length([minE:maxE]));
MLEP=zeros(length(w),length([minE:maxE]));
parfor ww=1:length(w)
    ptravel=ptravelB.*w(ww);
CPxTNS=ones(NS2,length([minE:maxE]));
PxTNS=ones(NS2,length([minE:maxE]));
    for ii=minE:maxE
       for mm=1:NS2
           f=find(TBNS(mm,:)>ii);
          g=find(E(mm,f)<=ii);
          dt=ptravel*(1-ptravel).^(ii-E(mm,f(g)));
          PxTNS(mm,ii-(minE)+1)=PxTNS(mm,ii-(minE)+1).*(1-prod(1-dt));
          gg=find(E(mm,:)<=ii);
          temps=min(TBNS(mm,gg),ii+1);
          dts=1-prod((1-ptravel).^temps);
          CPxTNS(mm,ii-(minE)+1)=CPxTNS(mm,ii-(minE)+1).*dts;
       end
    end
    MLE(ww,:)=mean(CPxTNS,1);
    MLEP(ww,:)=mean(PxTNS,1);
end

save('Weighted_Travel_Infectious_Country.mat','MLE','w','MLEP');