% Computes the time after arrival set a negative binomial with mean and
% std reportd in NEJM

% Serial interval (Time of symptom onset index patient to time of symptom
% onset in 
munb=7.5;
v=3.4^2;
p=munb/v;
r=munb.*p./(1-p);
NS1=10^4;
options=optimset('MaxIter',10^6,'TolFun',10^(-16),'TolX',10^(-16),'display','off');
[gp]=fmincon(@(x)((gaminv(0.025,x(1),5.2./x(1))-4.1).^2+(gaminv(0.975,x(1),5.2./x(1))-7).^2),100,[],[],[],[],0,1000,[],options);
gaminv(0.025,gp,5.2./gp)
gaminv(0.975,gp,5.2./gp)
mun=gamrnd(gp,5.2/gp,NS1,1);

Incu2=[0:21];
Test=zeros(22,1);
for ii=0:21
    t=0;
   for jj=ii:42
       t=t+(jj-ii).*nbinpdf(jj,r,p)./(1-nbincdf(ii-1,r,p));
   end
   Test(ii+1)=t;
end
[~,pdf] = IncubationDist(5.2,0);

MLE=pdf*Test;
UMLE=zeros(NS1,1);
for jj=1:NS1
    [~,pdf] = IncubationDist(mun(jj),0);
    UMLE(jj)=pdf*Test;
end

load('Probability_Travel.mat','F','pc');
w=exp(F)./sum(exp(F));
wc=cumsum(w);

ptravel=pc(F==max(F));

for ii=0:21
    t=0;
   for jj=0:(ii-1)
       t=t+(ii-jj).*ptravel.*(1-ptravel).^jj./(1-(1-ptravel).^ii);
   end
   Test(ii+1)=t;
end

[~,pdf] = IncubationDist(5.2,0);

MLET=pdf*Test;
UMLET=zeros(NS1,1);

for jj=1:NS1    
    [~,pdf] = IncubationDist(mun(jj),0);
    UMLET(jj)=pdf*Test;
end

NS2=10^6;
IND1=randi(NS1,NS2,1);
IND2=randi(NS1,NS2,1);
TtoT=UMLET(IND1)+UMLE(IND2);





