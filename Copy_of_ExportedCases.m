load('Par_NBin.mat','G','CV');
[Inc, IncO]=IncidenceData;
NS1=10^3;
NS2=10^2;
[R,Pr] = SampleDistShape(G,CV,NS1);
T=[];
for ii=1:length(Inc(:,1))
   T=[T Inc(ii,1).*ones(1,Inc(ii,2))]; 
end
for ii=1:length(IncO(:,1))
   T=[T IncO(ii,1).*ones(1,IncO(ii,2))]; 
end

INDX=datenum('01-22-2020')-datenum('12-29-2019')+1;
minE=-60;%min(E(:));
maxE=30;

TA=ones(1,length([minE:maxE]));
TAN=TA;
TA(INDX-minE:end)=0;
load('InfectiousPeriod.mat','P')
CDFP=cumsum(P(:,2))./sum(P(:,2)); % starts at zero

%% MLE Estimate

PI=zeros(NS2,length(T));
PIN=zeros(NS2,length(T));
PINS=zeros(NS2,length(T));

    IP=nbinrnd(G(1),G(2),NS2,length(T));
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
    TimeA(TimeA>30)=30;
    for ii=1:NS2
        for jj=1:length(T)
            if(E(ii,jj)<=(T(jj)-1)) % Subtract one as T indicates the time of symptoms. Thus for an incubation period of one day [E(ii,jj):(T(jj)-1)] needs to be length 1
                pt=0.005.*TA([E(ii,jj):(T(jj)-1)]-minE+1);
                PI(ii,jj)=LikelihoodMissed(pt);
                pt=0.005.*TAN([E(ii,jj):(T(jj)-1)]-minE+1);
                PIN(ii,jj)=LikelihoodMissed(pt);                
                pt=0.005.*TAN([E(ii,jj):(TimeA(ii,jj)-1)]-minE+1);
                PINS(ii,jj)=LikelihoodMissed(pt);
            end
        end
    end
MLExT=[min(T):max(T)];
MLExTN=[min(T):max(T)];
MLExTNS=[min(T):max(T)];
for ii=1:length(MLExT)
   f=find(T==MLExT(ii));
   if(~isempty(f))
        MLExT(ii)=sum(sum(PI(:,f)))./NS2;
        MLExTN(ii)=sum(sum(PIN(:,f)))./NS2;
        MLExTNS(ii)=sum(sum(PINS(:,f)))./NS2;
   else
       MLExT(ii)=0; 
       MLExTN(ii)=0;
       MLExTNS(ii)=0;
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
    TimeA(TimeA>30)=30;
    
    PI=zeros(NS2,length(T));
    PIN=zeros(NS2,length(T));
    PINS=zeros(NS2,length(T));
    for ii=1:NS2
        parfor jj=1:length(T)
            if(E(ii,jj)<=(T(jj)-1)) % Subtract one as T indicates the time of symptoms. Thus for an incubation period of one day [E(ii,jj):(T(jj)-1)] needs to be length 1
                pt=0.005.*TA([E(ii,jj):(T(jj)-1)]-minE+1);
                PI(ii,jj)=LikelihoodMissed(pt);
                pt=0.005.*TAN([E(ii,jj):(T(jj)-1)]-minE+1);
                PIN(ii,jj)=LikelihoodMissed(pt);                
                pt=0.005.*TAN([E(ii,jj):(TimeA(ii,jj)-1)]-minE+1);
                PINS(ii,jj)=LikelihoodMissed(pt);
            end
        end
    end
    parfor ii=1:length(UxT(nn,:))
       f=find(T==UxT(nn,ii));
       if(~isempty(f))
            UxT(nn,ii)=sum(sum(PI(:,f)))./NS2;
            UxTN(nn,ii)=sum(sum(PIN(:,f)))./NS2;
            UxTNS(nn,ii)=sum(sum(PINS(:,f)))./NS2;
       else
           UxT(nn,ii)=0; 
           UxTN(nn,ii)=0;
           UxTNS(nn,ii)=0;
       end
    end    
    UP(nn)=mean(sum(PI,2));
    UPN(nn)=mean(sum(PI,2));
    UPNS(nn)=mean(sum(PI,2));
end

save('Test_Travel_Ban.mat');