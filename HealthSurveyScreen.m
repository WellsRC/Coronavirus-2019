load('Weight_Flights.mat','FlightAll')
tf=strcmp({'China'},{FlightAll{:,1}});
wtc=1-[FlightAll{tf,2}];

NS1=10^4;
options=optimset('MaxIter',10^6,'TolFun',10^(-16),'TolX',10^(-16),'display','off');
[gp]=fmincon(@(x)((gaminv(0.025,x(1),5.2./x(1))-4.1).^2+(gaminv(0.975,x(1),5.2./x(1))-7).^2),100,[],[],[],[],0,1000,[],options);
gaminv(0.025,gp,5.2./gp)
gaminv(0.975,gp,5.2./gp)
mun=gamrnd(gp,5.2/gp,NS1,1);

Test=zeros(22,1);

load('Probability_Travel_Infection.mat','F','pc');
w=exp(F)./sum(exp(F));
wc=cumsum(w);

ptravel=pc(F==max(F));
Test(1)=1;

[~,pdf] = IncubationDist(5.2,0);
MLET=zeros(14,1);
for mm=1:14
    for ii=1:21
        t=0;
       for jj=0:(ii-1)
           if(jj<=mm)
            t=t+1.*wtc.*ptravel.*(1-ptravel).^jj./(wtc.*(1-(1-ptravel).^ii));
           end
       end
       Test(ii+1)=t;
    end
MLET(mm)=pdf*Test;
end


r=rand(NS1,1);
spc=zeros(NS1,1);

for ii=1:NS1
    f=find(r(ii)<=wc);
    f=f(1);
    spc(ii)=pc(f);
end

UMLET=zeros(NS1,14);
for ss=1:NS1    
    [~,pdf] = IncubationDist(mun(ss),0);
    for mm=1:14
        for ii=1:21
            t=0;
           for jj=0:(ii-1)
               if(jj<=mm)
                t=t+1.*wtc.*spc(ss).*(1-spc(ss)).^jj./(wtc.*(1-(1-spc(ss)).^ii));
               end
           end
           Test(ii+1)=t;
        end
    UMLET(ss,mm)=pdf*Test;
    end
end


save('Time_Screening.mat');





