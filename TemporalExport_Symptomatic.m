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


load('Probability_Travel.mat','F','pc');
w=exp(F)./sum(exp(F));
wc=cumsum(w);

ptravel=pc(F==max(F));


INDX=datenum('01-23-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
INDX2=datenum('01-25-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
INDXMV=datenum('01-1-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero

minE=-22;%min(E(:));
maxE=53;

TA=ones(1,length([minE:maxE]));
TAN=ones(1,length([minE:maxE]));
TA(INDX-minE:end)=0;

%% MLE Estimates

UxT=zeros(NS2,53);
IP=zeros(NS2,length(T)+length(TW)+length(TF));

% Sample incubation period
for ii=1:NS2
    [IP(ii,:),~] = IncubationDist(5.2,length(T)+length(TW)+length(TF));
end

% Travel Ban and screening
E=repmat([T TW TF],NS2,1)-IP; 
TT=[repmat([T],NS2,1) min(repmat([TW],NS2,1),INDX) min(repmat([TF],NS2,1),INDX2)];
D=(TT-E);
D(D<0)=0;

PI=1-(1-ptravel).^D;
for ii=1:53
   f=find([T TW TF]==ii);
   if(~isempty(f))
        UxT(:,ii)=sum(PI(:,f),2);
   end
end    

MLExTS=mean(UxT,1);
MLETS=mean(sum(PI,2)); 

% Screening and no travel ban
E=repmat([T TW TF],NS2,1)-IP; 
TT=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
D=(TT-E);
D(D<0)=0;

PI=1-(1-ptravel).^D;
for ii=1:53
   f=find([T TW TF]==ii);
   if(~isempty(f))
        UxT(:,ii)=sum(PI(:,f),2);
   end
end    

MLExS=mean(UxT,1);
MLES=mean(sum(PI,2)); 

% No screening and no travel ban

E=repmat([T TW TF],NS2,1)-IP; 
TT=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
TAO=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
TAOT=TAO;
TAO(TAOT>=INDXMV)=TT(TAOT>=INDXMV)+TimeMedJan1(length(TAOT(TAOT>=INDXMV))); 
TAO(TAOT<INDXMV)=TT(TAOT<INDXMV)+TimeMedDec31(length(TAOT(TAOT<INDXMV))); 
D=(TT-E);
D(D<0)=0;

Dt=(TAO-E);
Dt(Dt<0)=0;
PI=1-(1-ptravel).^Dt;

PItemp=1-(1-ptravel).^D;
TempT=[T TW TF];
for ii=1:53
   f=find(TempT==ii);
   if(~isempty(f))
        UxT(:,ii)=sum(PItemp(:,f),2);
   end
   for mm=1:NS2
      f=find(TAO(mm,:)>ii);
      g=find(TempT(f)<=ii);
      dt=ptravel*(1-ptravel).^(ii-TempT(f(g)));
      UxT(mm,ii)=UxT(mm,ii)+sum(dt);
   end
end    

MLEx=mean(UxT,1);
MLE=mean(sum(PI,2)); 


% No screening and travel ban

E=repmat([T TW TF],NS2,1)-IP; 
TT=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
TAO=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
TAOT=TAO;
TAO(TAOT>=INDXMV)=TT(TAOT>=INDXMV)+TimeMedJan1(length(TAOT(TAOT>=INDXMV))); 
TAO(TAOT<INDXMV)=TT(TAOT<INDXMV)+TimeMedDec31(length(TAOT(TAOT<INDXMV))); 
TBNS=TAO;
TBNS(:,(length(T)+1):(length(T)+length(TW)))=min(TBNS(:,(length(T)+1):(length(T)+length(TW))),INDX);
TBNS(:,(length(T)+length(TW)+1):end)=min(TBNS(:,(length(T)+length(TW)+1):end),INDX2);

D=(TT-E);
D(D<0)=0;

Dt=(TBNS-E);
Dt(Dt<0)=0;
PI=1-(1-ptravel).^Dt;

PItemp=1-(1-ptravel).^D;
TempT=[T TW TF];
for ii=1:53
   f=find(TempT==ii);
   if(~isempty(f))
        UxT(:,ii)=sum(PItemp(:,f),2);
   end
   for mm=1:NS2
      f=find(TBNS(mm,:)>ii);
      g=find(TempT(f)<=ii);
      dt=ptravel*(1-ptravel).^(ii-TempT(f(g)));
      UxT(mm,ii)=UxT(mm,ii)+sum(dt);
   end
end    

MLExNS=mean(UxT,1);
MLENS=mean(sum(PI,2)); 



%% Uncertainty

UMLExTS=zeros(NS1,length(MLExTS));
UMLExS=zeros(NS1,length(MLExTS));
UMLEx=zeros(NS1,length(MLExTS));
UMLExNS=zeros(NS1,length(MLExTS));


UMLETS=zeros(NS1,1);
UMLES=zeros(NS1,1);
UMLE=zeros(NS1,1);
UMLENS=zeros(NS1,1);
for ss=1:NS1
    UxT=zeros(NS2,53);
    IP=zeros(NS2,length(T)+length(TW)+length(TF));

    % Sample incubation period
    for ii=1:NS2
        [IP(ii,:),~] = IncubationDist(5.2,length(T)+length(TW)+length(TF));
    end

    % Travel Ban and screening
    E=repmat([T TW TF],NS2,1)-IP; 
    TT=[repmat([T],NS2,1) min(repmat([TW],NS2,1),INDX) min(repmat([TF],NS2,1),INDX2)];
    D=(TT-E);
    D(D<0)=0;

    PI=1-(1-ptravel).^D;
    for ii=1:53
       f=find([T TW TF]==ii);
       if(~isempty(f))
            UxT(:,ii)=sum(PI(:,f),2);
       end
    end    

    UMLExTS(ss,:)=mean(UxT,1);
    UMLETS(ss)=mean(sum(PI,2)); 

    % Screening and no travel ban
    E=repmat([T TW TF],NS2,1)-IP; 
    TT=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
    D=(TT-E);
    D(D<0)=0;

    PI=1-(1-ptravel).^D;
    for ii=1:53
       f=find([T TW TF]==ii);
       if(~isempty(f))
            UxT(:,ii)=sum(PI(:,f),2);
       end
    end    

    UMLExS(ss,:)=mean(UxT,1);
    UMLES(ss)=mean(sum(PI,2)); 

    % No screening and no travel ban

    E=repmat([T TW TF],NS2,1)-IP; 
    TT=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
    TAO=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
    TAOT=TAO;
    TAO(TAOT>=INDXMV)=TT(TAOT>=INDXMV)+TimeMedJan1(length(TAOT(TAOT>=INDXMV))); 
    TAO(TAOT<INDXMV)=TT(TAOT<INDXMV)+TimeMedDec31(length(TAOT(TAOT<INDXMV))); 
    D=(TT-E);
    D(D<0)=0;

    Dt=(TAO-E);
    Dt(Dt<0)=0;
    PI=1-(1-ptravel).^Dt;

    PItemp=1-(1-ptravel).^D;
    TempT=[T TW TF];
    for ii=1:53
       f=find(TempT==ii);
       if(~isempty(f))
            UxT(:,ii)=sum(PItemp(:,f),2);
       end
       for mm=1:NS2
          f=find(TAO(mm,:)>ii);
          g=find(TempT(f)<=ii);
          dt=ptravel*(1-ptravel).^(ii-TempT(f(g)));
          UxT(mm,ii)=UxT(mm,ii)+sum(dt);
       end
    end    

    UMLEx(ss,:)=mean(UxT,1);
    UMLE(ss)=mean(sum(PI,2)); 


    % No screening and travel ban

    E=repmat([T TW TF],NS2,1)-IP; 
    TT=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
    TAO=[repmat([T],NS2,1) repmat([TW],NS2,1) repmat([TF],NS2,1)];
    TAOT=TAO;
    TAO(TAOT>=INDXMV)=TT(TAOT>=INDXMV)+TimeMedJan1(length(TAOT(TAOT>=INDXMV))); 
    TAO(TAOT<INDXMV)=TT(TAOT<INDXMV)+TimeMedDec31(length(TAOT(TAOT<INDXMV))); 
    TBNS=TAO;
    TBNS(:,(length(T)+1):(length(T)+length(TW)))=min(TBNS(:,(length(T)+1):(length(T)+length(TW))),INDX);
    TBNS(:,(length(T)+length(TW)+1):end)=min(TBNS(:,(length(T)+length(TW)+1):end),INDX2);

    D=(TT-E);
    D(D<0)=0;

    Dt=(TBNS-E);
    Dt(Dt<0)=0;
    PI=1-(1-ptravel).^Dt;

    PItemp=1-(1-ptravel).^D;
    TempT=[T TW TF];
    for ii=1:53
       f=find(TempT==ii);
       if(~isempty(f))
            UxT(:,ii)=sum(PItemp(:,f),2);
       end
       for mm=1:NS2
          f=find(TBNS(mm,:)>ii);
          g=find(TempT(f)<=ii);
          dt=ptravel*(1-ptravel).^(ii-TempT(f(g)));
          UxT(mm,ii)=UxT(mm,ii)+sum(dt);
       end
    end    

    UMLExNS(ss,:)=mean(UxT,1);
    UMLENS(ss)=mean(sum(PI,2)); 


end
save('Test_Travel_Ban.mat');