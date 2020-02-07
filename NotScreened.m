%% Computes the general probabilityes of traveling in incubation or infectious period
clear;
% Determine the weight ofr China
load('Weight_Flights.mat','FlightAll')
tf=strcmp({'China'},{FlightAll{:,1}});
wtc=1-[FlightAll{tf,2}];
% Obtain the pdf for the incubation period
[~,Ipdf] = IncubationDist(5.2,0); 
% Obtain the probabaility of travel
load('Probability_Travel_Infection.mat','F','pc');
w=exp(F)./sum(exp(F)); % likelihood weight
wc=cumsum(w); % cumulative sum ofr sampling
ptravel=pc(F==max(F)); % Set travel probabaility
% Fit the gamma distribtuion for the sampling of the mean
options=optimset('MaxIter',10^6,'TolFun',10^(-16),'TolX',10^(-16),'display','off');
[gp]=fmincon(@(x)((gaminv(0.025,x(1),5.2./x(1))-4.1).^2+(gaminv(0.975,x(1),5.2./x(1))-7).^2),100,[],[],[],[],0,1000,[],options);
NS=10^4;
mun=gamrnd(gp,5.2/gp,NS,1);
% Get he pdf for the different distribution
UIpdf=zeros(NS,length(Ipdf));
for ii=1:NS
    [~,UIpdf(ii,:)] = IncubationDist(mun(ii),0);
end
% Obtaine the pdfs for the different measures 
load('TimetoMedVisituptoDec31.mat','D')
MP=D(:,2)./sum(D(:,2));
load('TimetoMedVisitJan1onward.mat','D')
MPA=D(:,2)./sum(D(:,2));
load('TimetoHospitaluptoDec31.mat','D')
HP=D(:,2)./sum(D(:,2));
load('TimetoHospitalJan1onward.mat','D')
HPA=D(:,2)./sum(D(:,2));

% Initialize 
DI=21; % suratino of incubation period max evalaution
L=zeros(NS,1);
LNS=zeros(NS,1);
LNSA=zeros(NS,1);
LNSH=zeros(NS,1);
LNSHA=zeros(NS,1);
LCT=zeros(NS,15);
LCTIP=zeros(NS,15);
tempL=zeros(DI+1,1);
tempNS=zeros(DI+1,DI+1);
tempNSA=zeros(DI+1,DI+1);
tempNSH=zeros(DI+1,DI+1);
tempNSHA=zeros(DI+1,DI+1);
MLECT=zeros(1,15);
MLECTIP=zeros(1,15);

% Calculate averages for mle
for s=0:DI    
   pt=ptravel.*ones(s,1);
   tempL(s+1)=Ipdf(s+1).*wtc.*LikelihoodMissed(pt);
   for ti=0:(length(MP)-1)
        pt=ptravel.*ones(s+ti,1);
        tempNS(s+1,ti+1)=Ipdf(s+1).*MP(ti+1).*wtc.*LikelihoodMissed(pt);
        tempNSA(s+1,ti+1)=Ipdf(s+1).*MPA(ti+1).*wtc.*LikelihoodMissed(pt);
   end
   for ti=0:(length(HP)-1)
        pt=ptravel.*ones(s+ti,1);
        tempNSH(s+1,ti+1)=Ipdf(s+1).*HP(ti+1).*wtc.*LikelihoodMissed(pt);
        tempNSHA(s+1,ti+1)=Ipdf(s+1).*HPA(ti+1).*wtc.*LikelihoodMissed(pt);
   end
   for cc=0:14
       pt=ptravel.*ones(min([s cc]),1);       
       MLECT(cc+1)=MLECT(cc+1)+Ipdf(s+1).*wtc.*LikelihoodMissed(pt); 
       for ti=0:(length(MPA)-1)
            pt=ptravel.*ones(min([s+ti cc]),1);
             MLECTIP(cc+1)=MLECTIP(cc+1)+Ipdf(s+1).*MPA(ti+1).*wtc.*LikelihoodMissed(pt);
       end
   end
end

MLE=sum(tempL);
MLENS=sum(tempNS(:));
MLENSA=sum(tempNSA(:));
MLENSH=sum(tempNSH(:));
MLENSHA=sum(tempNSHA(:));

% Sample travel probability
r=rand(NS,1);
spc=zeros(NS,1);

for ii=1:NS
    f=find(r(ii)<=wc);
    f=f(1);
    spc(ii)=pc(f);
end

% Run average for the uncertainty
for ii=1:NS
    IP=UIpdf(ii,:);
    for s=0:DI    
       pt=spc(ii).*ones(s,1);
       tempL(s+1)=IP(s+1).*wtc.*LikelihoodMissed(pt);       
       for ti=0:(length(MP)-1)
            pt=spc(ii).*ones(s+ti,1);
            tempNS(s+1,ti+1)=IP(s+1).*MP(ti+1).*wtc.*LikelihoodMissed(pt);
            tempNSA(s+1,ti+1)=IP(s+1).*MPA(ti+1).*wtc.*LikelihoodMissed(pt);
       end
       for ti=0:(length(HP)-1)
            pt=spc(ii).*ones(s+ti,1);
            tempNSH(s+1,ti+1)=IP(s+1).*HP(ti+1).*wtc.*LikelihoodMissed(pt);
            tempNSHA(s+1,ti+1)=IP(s+1).*HPA(ti+1).*wtc.*LikelihoodMissed(pt);
       end
        for cc=0:14
            pt=spc(ii).*ones(min([cc s]),1);
            LCT(ii,cc+1)=LCT(ii,cc+1)+IP(s+1).*wtc.*LikelihoodMissed(pt);  
               for ti=0:(length(MPA)-1)
                    pt=ptravel.*ones(min([s+ti cc]),1);
                     LCTIP(cc+1)=LCTIP(cc+1)+Ipdf(s+1).*MPA(ti+1).*wtc.*LikelihoodMissed(pt);
               end
        end
    end
    L(ii)=sum(tempL);
    LNS(ii)=sum(tempNS(:));
    LNSA(ii)=sum(tempNSA(:));
    LNSH(ii)=sum(tempNSH(:));
    LNSHA(ii)=sum(tempNSHA(:));
end
save('TravelDuringInfection.mat','L','LNS','LNSA','LNSH','LNSHA','LCT','MLE','MLENS','MLENSA','MLENSH','MLENSHA','MLECT','MLECTIP','LCTIP');

