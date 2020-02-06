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

load('Probability_Travel_Infection.mat','F','pc');
w=exp(F)./sum(exp(F));
wc=cumsum(w);

pmb=[pc(F==max(F)) 0 0];
r=rand(10^4,1);
spc=zeros(10^4,1);

for ii=1:10^4
    f=find(r(ii)<=wc);
    f=f(1);
    spc(ii)=pc(f);
end

pmb(2)=prctile(spc,2.5);
pmb(3)=prctile(spc,97.5);


F=zeros(3,length([1:maxE])); % only need to go from one to maxE as no infectious case is past one and we are looking at the symptomatic cases only
for mm=1:length(pmb)
    ptravel=pmb(mm);
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
    F(mm,:)=mean(UxT,1);
end
startDateofSim = datenum('12-06-2019');% Start date
XTL=datestr([startDateofSim+[0:4:(maxE-1)]],'mm-dd-yy');

figure('units','normalized','outerposition',[0 0 1 1]);

for ii=1:maxE
   patch(ii+[-0.35 0.35 0.35 -0.35], [F(2,ii) F(2,ii) F(3,ii) F(3,ii)],'k','LineStyle','none','Facealpha',0.3);hold on
   plot(ii+linspace(-0.35,0.35,2),F(1,ii).*ones(1,2),'k','LineWidth',2);    
   scatter(IncO(ii,1),IncO(ii,2),40,'r','filled');
end
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:4:maxE],'XTickLabel',XTL,'Xminortick','on','Yminortick','on');
xtickangle(45);
xlim([1 maxE+0.5])
ylim([0 30]);
ylabel({'Incidence'},'Fontsize',18);
xlabel({'Date'},'Fontsize',18);

plot([INDX INDX],[0 30],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);
plot([INDX2 INDX2],[0 30],'-.','color',[0.7 0.7 0.7],'LineWidth',1.5);