[R,Pr] = IncubationDist(5.2); % Avg 5.2 with CI of 4.1 to 7.0
[IncC,IncW,IncO]=IncidenceData;
NS1=10^3;
NS2=10^2;
T=[];
TW=[];
for ii=1:length(IncC(:,1))
   TW=[TW IncC(ii,1).*ones(1,IncC(ii,2))]; 
end
for ii=1:length(IncO(:,1))
   T=[T IncO(ii,1).*ones(1,IncO(ii,2))]; 
end

for ii=1:length(IncO(:,1))
   TW=[TW IncW(ii,1).*ones(1,IncW(ii,2))]; 
end

INDX=datenum('01-22-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
INDXMV=datenum('01-1-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
minE=-60;%min(E(:));
maxE=53;

TA=ones(1,length([minE:maxE]));
TAN=TA;
TA(INDX-minE:end)=0;


%% MLE Estimate
rng('shuffle');
PI=zeros(NS2,length(T)+length(TW));
PIN=zeros(NS2,length(T)+length(TW));
PINS=zeros(NS2,length(T)+length(TW));
PINSTr=zeros(NS2,length(T)+length(TW));

IP=nbinrnd(R,Pr,NS2,length(T)+length(TW));
E=repmat([T TW],NS2,1)-IP; 
TimeA=zeros(size(E));
for ii=1:NS2     
    TDEC30=TimeMedDec31(length(T)+length(TW));
    TJAN1=TimeMedDec31(length(T)+length(TW));
    for jj=1:length(T)
        if(T(jj)<INDXMV)
            TimeA(ii,jj)=T(jj)+TDEC30(jj);
        else
            TimeA(ii,jj)=T(jj)+TJAN1(jj);            
        end
    end
    for jj=(length(T)+1):(length(T)+length(TW))
        if(TW(jj-length(T))<INDXMV)
            TimeA(ii,jj)=TW(jj-length(T))+TDEC30(jj);
        else
            TimeA(ii,jj)=TW(jj-length(T))+TJAN1(jj);            
        end
    end
end
SimC=zeros(NS2,maxE);
for ii=1:NS2
    for jj=1:length(T)
        if(E(ii,jj)<=(T(jj)-1)) % Subtract one as T indicates the time of symptoms. Thus for an incubation period of one day [E(ii,jj):(T(jj)-1)] needs to be length 1
            pt=0.005.*TAN([E(ii,jj):(T(jj)-1)]-minE+1); % There is no lock-down for these affected areas
            PI(ii,jj)=LikelihoodMissed(pt);
            pt=0.005.*TAN([E(ii,jj):(T(jj)-1)]-minE+1);
            PIN(ii,jj)=LikelihoodMissed(pt);                     
            PINSTr(ii,jj)=LikelihoodMissed(pt); % Used in the temporal examination
            pt=0.005.*TAN([E(ii,jj):min([(TimeA(ii,jj)-1) maxE])]-minE+1);
            PINS(ii,jj)=LikelihoodMissed(pt);
        end
        if(T(jj)<TimeA(ii,jj))
            d=length([E(ii,jj):(T(jj)-1)]);
            pt=0.005.*TAN([T(jj):min([(TimeA(ii,jj)-1) maxE])]-minE+1);
            temp=(1-PINSTr(ii,jj)).*pt; % Probability not traveled during symptomatic period time probability of travel each day during syptoms
            for tt=2:length(temp)
               temp(tt)=temp(tt).*(1-LikelihoodMissed(pt(1:(tt-1))));
            end
            SimC(ii,T(jj):min([(TimeA(ii,jj)-1) maxE]))=SimC(ii,T(jj):min([(TimeA(ii,jj)-1) maxE]))+temp;
        end
    end
    for jj=(1+length(T)):(length(T)+length(TW))
        if(E(ii,jj)<=(TW(jj-length(T))-1)) % Subtract one as T indicates the time of symptoms. Thus for an incubation period of one day [E(ii,jj):(T(jj)-1)] needs to be length 1
            pt=0.005.*TA([E(ii,jj):(TW(jj-length(T))-1)]-minE+1); % lock down in effect
            PI(ii,jj)=LikelihoodMissed(pt);
            pt=0.005.*TAN([E(ii,jj):(TW(jj-length(T))-1)]-minE+1);
            PIN(ii,jj)=LikelihoodMissed(pt);                     
            PINSTr(ii,jj)=LikelihoodMissed(pt); % Used in the temporal examination
            pt=0.005.*TAN([E(ii,jj):min([(TimeA(ii,jj)-1) maxE])]-minE+1);
            PINS(ii,jj)=LikelihoodMissed(pt);
        end
        if(TW(jj-length(T))<TimeA(ii,jj))
            d=length([E(ii,jj):(TW(jj-length(T))-1)]);
            pt=0.005.*TAN([TW(jj-length(T)):min([(TimeA(ii,jj)-1) maxE])]-minE+1);
            temp=(1-PINSTr(ii,jj)).*pt; % Probability not traveled during symptomatic period time probability of travel each day during syptoms
            for tt=2:length(temp)
               temp(tt)=temp(tt).*(1-LikelihoodMissed(pt(1:(tt-1))));
            end
            SimC(ii,TW(jj-length(T)):min([(TimeA(ii,jj)-1) maxE]))=SimC(ii,TW(jj-length(T)):min([(TimeA(ii,jj)-1) maxE]))+temp;
        end
    end
end

MLExT=[1:maxE];
MLExTN=[1:maxE];
MLExTNS=[1:maxE];
for ii=1:length(MLExT)
   f=find([T TW]==MLExT(ii));
   if(~isempty(f))
        MLExT(ii)=sum(sum(PI(:,f),2))./NS2;
        MLExTN(ii)=sum(sum(PIN(:,f),2))./NS2;
        MLExTNS(ii)=sum(sum(PINSTr(:,f),2)+SimC(:,ii))./NS2;
   else
       MLExT(ii)=0; 
       MLExTN(ii)=0;
       MLExTNS(ii)=sum(SimC(:,ii))./NS2;
   end
end
MLEP=mean(sum(PI,2));
MLEPN=mean(sum(PIN,2));
MLEPNS=mean(sum(PINS,2));

%% Uncertainty

UxT=repmat([min(T):max(T)],NS1,1);
UxTN=repmat([min(T):max(T)],NS1,1);
UxTNS=repmat([min(T):max(T)],NS1,1);
UP=zeros(NS1,1);
UPN=zeros(NS1,1);
UPNS=zeros(NS1,1);
for nn=1:NS1    
    G1=R(nn);
    G2=Pr(nn);
    IP=nbinrnd(G1,G2,NS2,length(T));
    E=repmat(T,NS2,1)-IP; 
    TimeA=zeros(size(E));
    rt=rand(NS2,length(T));
    for ii=1:NS2 
        for jj=1:length(T)
           ff=find(rt(ii,jj)<=CDFP);
           ff=ff(end)-1;
           TimeA(ii,jj)=T(jj)+ff;
        end
    end
    
    PI=zeros(NS2,length(T));
    PIN=zeros(NS2,length(T));
    PINS=zeros(NS2,length(T));
    PINSTr=zeros(NS2,length(T));
    SimC=zeros(NS2,length([min(T):max(T)]));
    
    for ii=1:NS2
        for jj=1:length(T)
            if(E(ii,jj)<=(T(jj)-1)) % Subtract one as T indicates the time of symptoms. Thus for an incubation period of one day [E(ii,jj):(T(jj)-1)] needs to be length 1
                pt=0.005.*TA([E(ii,jj):(T(jj)-1)]-minE+1);
                PI(ii,jj)=LikelihoodMissed(pt);
                pt=0.005.*TAN([E(ii,jj):(T(jj)-1)]-minE+1);
                PIN(ii,jj)=LikelihoodMissed(pt);      
                PINSTr(ii,jj)=LikelihoodMissed(pt); % Used in the temporal examination
                pt=0.005.*TAN([E(ii,jj):min([(TimeA(ii,jj)-1) maxE])]-minE+1);
                PINS(ii,jj)=LikelihoodMissed(pt);
            end
            
            if(T(jj)<TimeA(ii,jj))
                d=length([E(ii,jj):(T(jj)-1)]);
                pt=0.005.*TAN([T(jj):min([(TimeA(ii,jj)-1) maxE])]-minE+1);
                temp=(1-PINSTr(ii,jj)).*pt; % Probability not traveled during symptomatic period time probability of travel each day during syptoms
                for tt=2:length(temp)
                   temp(tt)=temp(tt).*(1-LikelihoodMissed(pt(1:(tt-1))));
                end
                SimC(ii,T(jj):min([(TimeA(ii,jj)-1) maxE]))=SimC(ii,T(jj):min([(TimeA(ii,jj)-1) maxE]))+temp;
            end
        end
    end
    parfor ii=1:length(UxT(nn,:))
       f=find(T==UxT(nn,ii));
       if(~isempty(f))
            UxT(nn,ii)=sum(sum(PI(:,f)))./NS2;
            UxTN(nn,ii)=sum(sum(PIN(:,f)))./NS2;
            UxTNS(nn,ii)=sum(sum(PINSTr(:,f),2)+SimC(:,ii))./NS2;
       else
           UxT(nn,ii)=0; 
           UxTN(nn,ii)=0;
           UxTNS(nn,ii)=sum(SimC(:,ii))./NS2;
       end
    end    
    UP(nn)=mean(sum(PI,2));
    UPN(nn)=mean(sum(PIN,2));
    UPNS(nn)=mean(sum(PINS,2));
end

save('Test_Travel_Ban.mat');