%% Use of the serial interval for non-screening
[~,Ipdf] = IncubationDist(5.2,0);
load('Probability_Travel.mat','F','pc');
w=exp(F)./sum(exp(F));
wc=cumsum(w);

ptravel=pc(F==max(F));
options=optimset('MaxIter',10^6,'TolFun',10^(-16),'TolX',10^(-16),'display','off');
[gp]=fmincon(@(x)((gaminv(0.025,x(1),5.2./x(1))-4.1).^2+(gaminv(0.975,x(1),5.2./x(1))-7).^2),100,[],[],[],[],0,1000,[],options);
gaminv(0.025,gp,5.2./gp)
gaminv(0.975,gp,5.2./gp)
NS=10^4;
mun=gamrnd(gp,5.2/gp,NS,1);
UIpdf=zeros(NS,length(Ipdf));
for ii=1:NS
    [~,UIpdf(ii,:)] = IncubationDist(mun(ii),0);
end
load('TimetoMedVisituptoDec31.mat','D')
MP=D(:,2)./sum(D(:,2));
load('TimetoMedVisitJan1onward.mat','D')
MPA=D(:,2)./sum(D(:,2));
load('TimetoHospitaluptoDec31.mat','D')
HP=D(:,2)./sum(D(:,2));
load('TimetoHospitalJan1onward.mat','D')
HPA=D(:,2)./sum(D(:,2));
DI=21;
tempL=zeros(DI+1,1);
w=[0:0.005:1];
MLE=zeros(length(w),1);
for ii=1:length(w)
    L=zeros(NS,1);
    for s=0:DI    
       pt=ptravel.*ones(s,1).*w(ii);
       tempL(s+1)=Ipdf(s+1).*LikelihoodMissed(pt);
    end
    MLE(ii)=sum(tempL);
end

save('Weighted_Travel_Inubation.mat','MLE','w');