load('Par_NBin.mat','G','GU');
load('Incidence_Start_Dec292019.mat')
load('CumulativeIncidence.mat','Other')
INDX=datenum('01-21-2020')-datenum('12-29-2019')+1;
TempI=[China(2:end)-China(1:end-1)]'; % Starts Jan. 21 as index starts the 20th
TempI=[[INDX:(INDX+length(TempI)-1)]' TempI];
for ii=1:length(Inc(:,1))
    ff=find(Inc2(:,1)==Inc(ii,1));
    if(~isempty(ff))
        Inc(ii,2)=max([Inc(ii,2) Inc2(ff,2)]);
    end
    ff=find(TempI(:,1)==Inc(ii,1));
    if(~isempty(ff))
        Inc(ii,2)=max([Inc(ii,2) TempI(ff,2)]);
    end
end
NS=10^3;
E=zeros(NS,sum(Inc(:,2))+sum(IncO(:,2)));
T=[];
for ii=1:length(Inc(:,1))
   T=[T Inc(ii,1).*ones(1,Inc(ii,2))]; 
end
for ii=1:length(IncO(:,1))
   T=[T IncO(ii,1).*ones(1,IncO(ii,2))]; 
end

INDX=datenum('01-22-2020')-datenum('12-29-2019')+1;
minE=-60;
maxE=max(E(:));
RR=zeros(NS,length(minE:maxE));
RRTB=zeros(NS,length(minE:maxE));

TA=ones(1,length([minE:maxE]));
TAN=TA;
TA(INDX-minE:end)=0;

for ii=1:NS
    for jj=1:length(T)
       RR(ii,[[E(ii,jj):T(jj)]-minE+1])= RR(ii,[[E(ii,jj):T(jj)]-minE+1])+1;
       RRTB(ii,[[E(ii,jj):T(jj)]-minE+1])= RRTB(ii,[[E(ii,jj):T(jj)]-minE+1])+TA([[E(ii,jj):T(jj)]-minE+1]);
    end
end

PR=1-(1-0.005).^RR;
PRTB=1-(1-0.005).^RRTB;
plot([minE:maxE],mean(PR,1),'r',[minE:maxE],mean(PRTB,1),'k')