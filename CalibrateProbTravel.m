clear;
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
options=optimset('MaxIter',10^6,'TolFun',10^(-16),'TolX',10^(-16),'display','off');
[gp]=fmincon(@(x)((gaminv(0.025,x(1),5.2./x(1))-4.1).^2+(gaminv(0.975,x(1),5.2./x(1))-7).^2),100,[],[],[],[],0,1000,[],options);
gaminv(0.025,gp,5.2./gp)
gaminv(0.975,gp,5.2./gp)
mun=gamrnd(gp,5.2/gp,NS1,1);


INDX=datenum('01-23-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
INDX2=datenum('01-25-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero
INDXMV=datenum('01-1-2020')-datenum('12-06-2019')+1; % Need to add one since the week index for Dec 6 would be zero

minE=-22;%min(E(:));
maxE=53;



pc=0.001.*[1.3:0.01:5];%linspace(1.3,3.5,101);
F=zeros(length(pc),1);
for mm=1:length(pc)
    UxT=zeros(NS1*NS2,53);
    IP=zeros(NS1*NS2,length(T)+length(TW)+length(TF));
    for nn=1:NS1    
        for ii=1:NS2
            [IP(ii+NS2.*(nn-1),:),~] = IncubationDist(5.2,length(T)+length(TW)+length(TF));
        end
    end
    E=repmat([T TW TF],NS2*NS1,1)-IP; 
    TT=[repmat([T],NS2*NS1,1)-1 min(repmat([TW],NS2*NS1,1)-1,INDX-1) min(repmat([TF],NS2*NS1,1)-1,INDX2-1)];
    D=(TT-E)+1;
    D(D<0)=0;

    PI=1-(1-pc(mm)).^D;
    for ii=1:53
       f=find([T TW TF]==ii);
       if(~isempty(f))
            UxT(:,ii)=sum(PI(:,f),2);
       end
    end    
    M=poisspdf(IncOutside',mean(UxT,1));
    F(mm)=sum(log(M));
end
save('Probability_Travel.mat','F','pc');

plot(pc,F);