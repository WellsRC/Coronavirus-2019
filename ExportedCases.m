% load('Par_NBin.mat','G','GU');
% load('Incidence_Start_Dec292019.mat')
% load('CumulativeIncidence.mat')
% INDX=datenum('01-21-2020')-datenum('12-29-2019')+1;
% TempI=[China(2:end)-China(1:end-1)]'; % Starts Jan. 21 as index starts the 20th
% TempI=[[INDX:(INDX+length(TempI)-1)]' TempI];
% for ii=1:length(Inc(:,1))
%     ff=find(Inc2(:,1)==Inc(ii,1));
%     if(~isempty(ff))
%         Inc(ii,2)=max([Inc(ii,2) Inc2(ff,2)]);
%     end
%     ff=find(TempI(:,1)==Inc(ii,1));
%     if(~isempty(ff))
%         Inc(ii,2)=max([Inc(ii,2) TempI(ff,2)]);
%     end
% end
% Inc=[Inc;TempI(end-3:end,:)];
% TempI=[Other(2:end)-Other(1:end-1)]'; % Starts Jan. 21 as index starts the 20th
% TempI=[[INDX:(INDX+length(TempI)-1)]' TempI];
% 
% for ii=1:length(IncO(:,1))
%     ff=find(TempI(:,1)==IncO(ii,1));
%     if(~isempty(ff))
%         IncO(ii,2)=max([IncO(ii,2) TempI(ff,2)]);
%     end
% end
% IncO=[IncO; TempI(end-1:end,:)];
% NS1=10^3;
% NS2=10^2;
% T=[];
% for ii=1:length(Inc(:,1))
%    T=[T Inc(ii,1).*ones(1,Inc(ii,2))]; 
% end
% for ii=1:length(IncO(:,1))
%    T=[T IncO(ii,1).*ones(1,IncO(ii,2))]; 
% end
% 
% INDX=datenum('01-22-2020')-datenum('12-29-2019')+1;
% minE=-100;%min(E(:));
% maxE=30;
% 
% TA=ones(1,length([minE:maxE]));
% TAN=TA;
% TA(INDX-minE:end)=0;
% load('InfectiousPeriod.mat','P')
% CDFP=cumsum(P(:,2))./sum(P(:,2)); % starts at zero
% 
% PI=zeros(NS1*NS2,length(T));
% PIN=zeros(NS2*NS1,length(T));
% PINS=zeros(NS2*NS1,length(T));
% for nn=1:NS1    
%     G1=GU(1,1)+(GU(2,1)-GU(1,1)).*rand(1);
%     G2=GU(1,2)+(GU(2,2)-GU(1,2)).*rand(1);
%     IP=nbinrnd(G1,G2,NS2,length(T));
%     E=repmat(T,NS2,1)-IP; 
%     TimeA=zeros(size(E));
%     rt=rand(NS2,length(T));
%     for ii=1:NS2 
%         for jj=1:length(T)
%            ff=find(rt(ii,jj)<=CDFP);
%            ff=ff(end)-1;
%            TimeA(ii,jj)=T(jj)+ff;
%         end
%     end
%     TimeA(TimeA>30)=30;
%     for ii=1:NS2
%         for jj=1:length(T)
%             if(E(ii,jj)<=(T(jj)-1)) % Subtract one as T indicates the time of symptoms. Thus for an incubation period of one day [E(ii,jj):(T(jj)-1)] needs to be length 1
%                 pt=0.005.*TA([E(ii,jj):(T(jj)-1)]-minE+1);
%                 PI(ii+NS2.*(nn-1),jj)=LikelihoodMissed(pt);
%                 pt=0.005.*TAN([E(ii,jj):(T(jj)-1)]-minE+1);
%                 PIN(ii+NS2.*(nn-1),jj)=LikelihoodMissed(pt);                
%                 pt=0.005.*TAN([E(ii,jj):(TimeA(ii,jj)-1)]-minE+1);
%                 PINS(ii+NS2.*(nn-1),jj)=LikelihoodMissed(pt);
%             end
%         end
%     end
% end
% close all;
% MP=sum(PI,2);
% MPN=sum(PIN,2);
% MPNS=sum(PINS,2);
% save('Test_Travel_Ban.mat');
load('Test_Travel_Ban.mat');
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.052794117647059,0.632469428147312,0.283865546218487,0.341162790697676]); 
histogram(MP,[0:5:450],'Facecolor','b','LineStyle','none','FaceAlpha',0.3,'Normalization','probability'); hold on
histogram(MPN,[0:5:450],'Facecolor','k','LineStyle','none','FaceAlpha',0.3,'Normalization','probability'); hold on
histogram(MPNS,[0:5:450],'Facecolor','r','LineStyle','none','FaceAlpha',0.3,'Normalization','probability'); hold on
xlabel('Expected number of cases','Fontsize',18);
yh=ylabel('Frequency','Fontsize',18);
legend('Travel restriction','No travel restriction','No Restrictions and no screening');
legend boxoff;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16);
box off;