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


minE=1;%min(E(:));
maxE=53;

%% MLE Estimates
IP=zeros(NS2,length(T)+length(TW)+length(TF));

% Sample incubation period
for ii=1:NS2
    [IP(ii,:),~] = IncubationDist(5.2,length(T)+length(TW)+length(TF));
end

% Travel Ban and screening
E=repmat([T TW TF],NS2,1)-IP; 
TNR=repmat([T TW TF],NS2,1);

w=[0:0.005:1];
MLE=zeros(length(w),length([minE:maxE]));
for ww=1:length(w)
    
CPxS=ones(NS2,length([minE:maxE]));
    for ii=minE:maxE
       for mm=1:NS2
          gg=find(E(mm,:)<=ii);
          temps=min(TNR(mm,gg),ii+1);
          dts=1-prod((1-ptravel*(w(ww))).^temps);
          CPxS(mm,ii-(minE)+1)=CPxS(mm,ii-(minE)+1).*dts;
       end
    end
    MLE(ww,:)=mean(CPxS,1);
end

save('Weighted_Travel_Inubation.mat','MLE','w');