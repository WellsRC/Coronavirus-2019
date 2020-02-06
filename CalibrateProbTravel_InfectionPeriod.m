clear;
load('Weight_Flights.mat','FlightAll')
tf=strcmp({'China'},{FlightAll{:,1}});
wtc=1-[FlightAll{tf,2}]; % not to china
[IncC,IncW,IncH,IncO]=IncidenceData;
NS1=1;
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


INDX=datenum('01-23-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
INDX2=datenum('01-25-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
INDXMV=datenum('01-1-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero

minE=-22;%min(E(:));
maxE=max(IncO(:,1));


IP=zeros(NS1*NS2,length(T)+length(TW)+length(TF));
for nn=1:NS1    
    for ii=1:NS2
        [IP(ii+NS2.*(nn-1),:),~] = IncubationDist(5.2,length(T)+length(TW)+length(TF));
    end
end
E=repmat([T TW TF],NS2,1)-IP; 
TT=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
TAO=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
TAOT=TAO;
TAO(TAOT>=INDXMV)=TT(TAOT>=INDXMV)+TimeMedJan1(length(TAOT(TAOT>=INDXMV))); 
TAO(TAOT<INDXMV)=TT(TAOT<INDXMV)+TimeMedDec31(length(TAOT(TAOT<INDXMV))); 
TBNS=TAO;
TBNS(:,(length(T)+1):(length(T)+length(TW)))=min(TBNS(:,(length(T)+1):(length(T)+length(TW))),INDX);
TBNS(:,(length(T)+length(TW)+1):end)=min(TBNS(:,(length(T)+length(TW)+1):end),INDX2);

pc=0.01.*[0.3:0.01:3];%linspace(1.3,3.5,101);
F=zeros(length(pc),1);
parfor mm=1:length(pc)
    ptravel=pc(mm);
    UxT=zeros(NS1*NS2,maxE);

    D=(TT-E);
    D(D<0)=0;

    PItemp=wtc.*(1-(1-ptravel).^D);
    
    TempT=[T TW TF];
    for ii=1:maxE
       f=find(TempT==ii);
       if(~isempty(f))
            UxT(:,ii)=sum(PItemp(:,f),2);
       end
       for zz=1:NS2
          f=find(TBNS(zz,:)>ii);
          g=find(TempT(f)<=ii);
          dt=wtc.*ptravel*(1-ptravel).^(ii-TempT(f(g)));
          UxT(zz,ii)=UxT(zz,ii)+sum(dt);
       end
    end    
    M=poisspdf(IncOutside',mean(UxT,1));
    F(mm)=sum(log(M));
end
plot(pc,F)
save('Probability_Travel_Infection.mat','F','pc');
