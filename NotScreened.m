%% Use of the serial interval for non-screening
load('Par_NBin.mat','G','CV');

NS=10^4;
[R,Pr] = SampleDistShape(G,CV,NS);
load('Par_NBin_Serial.mat','M','CVM');

[RM,PrM] = SampleDistShape(M,CVM,NS);
DI=42;
L=zeros(NS,1);
LNS=zeros(NS,1);
LCT=zeros(NS,36);
LNSCT=zeros(NS,36);
tempL=zeros(DI+1,1);
tempNS=zeros(DI+1,1);
MLECT=zeros(1,36);
MLECTNS=zeros(1,36);
for s=0:DI    
   pt=0.005.*ones(s,1);
   tempL(s+1)=nbinpdf(s,G(1),G(2)).*LikelihoodMissed(pt);
   tempNS(s+1)=nbinpdf(s,M(1),M(2)).*LikelihoodMissed(pt);
   for cc=0:35
       pt=0.005.*ones(min([s cc]),1);       
       MLECT(cc+1)=MLECT(cc+1)+nbinpdf(s,G(1),G(2)).*LikelihoodMissed(pt);
       MLECTNS(cc+1)=MLECTNS(cc+1)+nbinpdf(s,M(1),M(2)).*LikelihoodMissed(pt);
   end
end

MLE=sum(tempL);
MLENS=sum(tempNS(:));


for ii=1:NS
    tempL=zeros(DI+1,1);
    G1=R(ii);
    G2=Pr(ii);
    M1=RM(ii);
    M2=PrM(ii);
    for s=0:DI    
       pt=0.005.*ones(s,1);
       tempL(s+1)=nbinpdf(s,G1,G2).*LikelihoodMissed(pt);
       tempNS(s+1)=nbinpdf(s,M1,M2).*LikelihoodMissed(pt);
        for cc=0:35
            pt=0.005.*ones(min([cc s]),1);
            LCT(ii,cc+1)=LCT(ii,cc+1)+nbinpdf(s,G1,G2).*LikelihoodMissed(pt);
            LNSCT(ii,cc+1)=LNSCT(ii,cc+1)+nbinpdf(s,M1,M2).*LikelihoodMissed(pt);
        end
    end
    L(ii)=sum(tempL);
    LNS(ii)=sum(tempNS(:));
end
save('Missed_Screening.mat','L','LNS','LCT','LNSCT','MLE','MLENS','MLECT','MLECTNS');
% close all;
% load('Missed_Screening.mat','L','LNS','LCT','LNSCT','MLE','MLENS','MLECT','MLECTNS');
% startDateofSim = datenum('12-29-2020');% Start date
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% subplot('Position',[0.052794117647059,0.632469428147312,0.283865546218487,0.341162790697676]); 
% histogram(LNS,[0:0.001:0.11],'Facecolor','r','LineStyle','none','FaceAlpha',0.3,'Normalization','probability'); hold on
% histogram(L,[0:0.001:0.11],'Facecolor','k','LineStyle','none','FaceAlpha',0.3,'Normalization','probability'); hold on
% xlabel('Probability infected person missed','Fontsize',18);
% yh=ylabel('Frequency','Fontsize',18);
% legend('No screening','Screening');
% legend boxoff;
% xlim([0 0.11]);
% ylim([0 0.04]);
% box off;
% set(gca,'LineWidth',2,'tickdir','out','Fontsize',16);
% text(yh.Extent(1),max(ylim),'A','Fontsize',32,'FontWeight','bold');
% [Inc, IncO]=IncidenceData;
% Inc(:,2)=cumsum(Inc(:,2));
% IncO(:,2)=cumsum(IncO(:,2));
% XTL=datestr([startDateofSim+[0:(length(Inc)-1)]],'mmm-dd');
% subplot('Position',[0.397666849348435,0.630443085695437,0.588677688466691,0.341162790697676]); 
% for ii=1:length(Inc(:,1))
%    patch(ii+[-0.35 -0.35 0.35 0.35],prctile(L,[2.5 97.5 97.5 2.5]).*(Inc(ii,2)+IncO(ii,2)),'k','Facealpha',0.3,'LineStyle','none'); hold on
%    plot(ii+linspace(-0.35,0.35,101),MLE.*(Inc(ii,2)+IncO(ii,2)).*ones(101,1),'k','LineWidth',2); hold on
%    
%    
%    patch(ii+[-0.35 -0.35 0.35 0.35],prctile(LNS,[2.5 97.5 97.5 2.5]).*(Inc(ii,2)+IncO(ii,2)),'r','Facealpha',0.3,'LineStyle','none'); hold on
%    plot(ii+linspace(-0.35,0.35,101),MLENS.*(Inc(ii,2)+IncO(ii,2)).*ones(101,1),'r','LineWidth',2); hold on
% end
% xlim([0.5 length(Inc(:,1))+0.5]);
% scatter([1:length(IncO(:,2))],IncO(:,2),40,'k','filled');
% yh=ylabel({'Expected number of','cases missed'},'Fontsize',18);
% xlabel('Date of epidemic (Year 2020)','Fontsize',18);
% box off;
% set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:length(Inc)],'XTickLabel',XTL,'YTick',[0:50:400]);
% xtickangle(45);
% ylim([0 400]);
% text(yh.Extent(1),max(ylim),'B','Fontsize',32,'FontWeight','bold');
% 
% subplot('Position',[0.052794117647059,0.17,0.283865546218487,0.341162790697676]); 
% patch([0:35 35:-1:0],[prctile(LCT,2.5) flip(prctile(LCT,97.5))],'k','Facealpha',0.3,'LineStyle','none'); hold on
% plot([0:35],MLECT,'k','LineWidth',2);
% 
% patch([0:35 35:-1:0],[prctile(LNSCT,2.5) flip(prctile(LNSCT,97.5))],'r','Facealpha',0.3,'LineStyle','none'); hold on
% plot([0:35],MLECTNS,'r','LineWidth',2);
% 
% xlim([0 35]);
% xlabel('Days from exposure to isolation','Fontsize',18)
% yh=ylabel('Expected probability','Fontsize',18);
% box off;
% set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[0:6:42],'Xminortick','on');
% text(yh.Extent(1),max(ylim),'C','Fontsize',32,'FontWeight','bold');
% 
% [Inc, IncO]=IncidenceData;
% 
% Inc(:,2)=cumsum(Inc(:,2));
% IncO(:,2)=cumsum(IncO(:,2));
% subplot('Position',[0.397666849348435,0.17,0.588677688466691,0.341162790697676]); 
% for ii=1:4
%    patch(ii+[-0.35 -0.35 0.35 0.35],[prctile(LCT(:,[3.*ii+1]),[2.5 97.5 97.5 2.5])].*(Inc(end,2)+IncO(end,2)),'k','Facealpha',0.3,'LineStyle','none'); hold on
%    plot(ii+linspace(-0.35,0.35,101),MLECT(3.*ii+1).*(Inc(end,2)+IncO(end,2)).*ones(101,1),'k','LineWidth',2); hold on   
%    
%    patch(ii+[-0.35 -0.35 0.35 0.35],[prctile(LNSCT(:,[3.*ii+1]),[2.5 97.5 97.5 2.5])].*(Inc(end,2)+IncO(end,2)),'r','Facealpha',0.3,'LineStyle','none'); hold on
%    plot(ii+linspace(-0.35,0.35,101),MLECTNS(3.*ii+1).*(Inc(end,2)+IncO(end,2)).*ones(101,1),'r','LineWidth',2); hold on
% end
% xh=xlabel('Days from exposure to isolation','Fontsize',18);
% xh.Position=[2.500001907348632,-37.18,-1];
% yh=ylabel({'Expected number of','cases missed'},'Fontsize',18);
% box off;
% 
% ylim([0 300]);
% xlim([0.5 4.5]);
% set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1 2 3 4],'XTickLabel',{'3','6','9','12'},'YTick',[0:25:300]);
% xtickangle(45);
% text(yh.Extent(1),max(ylim),'D','Fontsize',32,'FontWeight','bold');


% Use of the admission time to hospital from symptom onset for non-screening


load('InfectiousPeriod.mat','P')
MP=P(:,2)./sum(P(:,2));
NS=10^4;
DI=42;
L=zeros(NS,1);
LNS=zeros(NS,1);
LCT=zeros(NS,36);
LNSCT=zeros(NS,36);
tempL=zeros(DI+1,1);
tempNS=zeros(DI+1,DI+1);
MLECT=zeros(1,36);
MLECTNS=zeros(1,36);
for s=0:DI    
   pt=0.005.*ones(s,1);
   tempL(s+1)=nbinpdf(s,G(1),G(2)).*LikelihoodMissed(pt);
   for ti=0:(length(MP)-1)
        pt=0.005.*ones(s+ti,1);
        tempNS(s+1,ti+1)=nbinpdf(s,G(1),G(2)).*MP(ti+1).*LikelihoodMissed(pt);
   end
   for cc=0:35
       pt=0.005.*ones(min([s cc]),1);       
       MLECT(cc+1)=MLECT(cc+1)+nbinpdf(s,G(1),G(2)).*LikelihoodMissed(pt);       
       for ti=0:(length(MP)-1)
           pt=0.005.*ones(min([s+ti cc]),1);      
           MLECTNS(cc+1)=MLECTNS(cc+1)+nbinpdf(s,G(1),G(2)).*MP(ti+1).*LikelihoodMissed(pt);
       end
   end
end

MLE=sum(tempL);
MLENS=sum(tempNS(:));


for ii=1:NS
    tempL=zeros(DI+1,1);
    G1=R(ii);
    G2=Pr(ii);
    for s=0:DI    
       pt=0.005.*ones(s,1);
       tempL(s+1)=nbinpdf(s,G1,G2).*LikelihoodMissed(pt);       
       for ti=0:(length(MP)-1)
            pt=0.005.*ones(s+ti,1);
            tempNS(s+1,ti+1)=nbinpdf(s,G1,G2).*MP(ti+1).*LikelihoodMissed(pt);
       end
        for cc=0:35
            pt=0.005.*ones(min([cc s]),1);
            LCT(ii,cc+1)=LCT(ii,cc+1)+nbinpdf(s,G1,G2).*LikelihoodMissed(pt); 
           for ti=0:(length(MP)-1)
               pt=0.005.*ones(min([s+ti cc]),1);  
               LNSCT(ii,cc+1)=LNSCT(ii,cc+1)+nbinpdf(s,G1,G2).*MP(ti+1).*LikelihoodMissed(pt);
           end
        end
    end
    L(ii)=sum(tempL);
    LNS(ii)=sum(tempNS(:));
end
save('Missed_Screening_TimetoIsolation.mat','L','LNS','LCT','LNSCT','MLE','MLENS','MLECT','MLECTNS');
% close all;
% load('Missed_Screening_TimetoIsolation.mat','L','LNS','LCT','LNSCT','MLE','MLENS','MLECT','MLECTNS');
% startDateofSim = datenum('12-29-2020');% Start date
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% subplot('Position',[0.052794117647059,0.632469428147312,0.283865546218487,0.341162790697676]); 
% histogram(LNS,[0:0.001:0.11],'Facecolor','r','LineStyle','none','FaceAlpha',0.3,'Normalization','probability'); hold on
% histogram(L,[0:0.001:0.11],'Facecolor','k','LineStyle','none','FaceAlpha',0.3,'Normalization','probability'); hold on
% xlabel('Probability infected person missed','Fontsize',18);
% yh=ylabel('Frequency','Fontsize',18);
% legend('No screening','Screening');
% legend boxoff;
% xlim([0 0.11]);
% ylim([0 0.04]);
% box off;
% set(gca,'LineWidth',2,'tickdir','out','Fontsize',16);
% text(yh.Extent(1),max(ylim),'A','Fontsize',32,'FontWeight','bold');
% [Inc, IncO]=IncidenceData;
% 
% Inc(:,2)=cumsum(Inc(:,2));
% IncO(:,2)=cumsum(IncO(:,2));
% XTL=datestr([startDateofSim+[0:(length(Inc(:,2))-1)]],'mmm-dd');
% subplot('Position',[0.397666849348435,0.630443085695437,0.588677688466691,0.341162790697676]); 
% for ii=1:length(Inc(:,2))
%    patch(ii+[-0.35 -0.35 0.35 0.35],prctile(L,[2.5 97.5 97.5 2.5]).*(Inc(ii,2)+IncO(ii,2)),'k','Facealpha',0.3,'LineStyle','none'); hold on
%    plot(ii+linspace(-0.35,0.35,101),MLE.*(Inc(ii,2)+IncO(ii,2)).*ones(101,1),'k','LineWidth',2); hold on
%    
%    
%    patch(ii+[-0.35 -0.35 0.35 0.35],prctile(LNS,[2.5 97.5 97.5 2.5]).*(Inc(ii,2)+IncO(ii,2)),'r','Facealpha',0.3,'LineStyle','none'); hold on
%    plot(ii+linspace(-0.35,0.35,101),MLENS.*(Inc(ii,2)+IncO(ii,2)).*ones(101,1),'r','LineWidth',2); hold on
% end
% xlim([0.5 length(Inc(:,2))+0.5]);
% scatter([1:length(IncO(:,2))],IncO(:,2),40,'k','filled');
% yh=ylabel({'Expected number of','cases missed'},'Fontsize',18);
% xlabel('Date of epidemic (Year 2020)','Fontsize',18);
% box off;
% set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:length(Inc(:,2))],'XTickLabel',XTL,'YTick',[0:50:400]);
% xtickangle(45);
% ylim([0 400]);
% text(yh.Extent(1),max(ylim),'B','Fontsize',32,'FontWeight','bold');
% 
% subplot('Position',[0.052794117647059,0.17,0.283865546218487,0.341162790697676]); 
% patch([0:35 35:-1:0],[prctile(LCT,2.5) flip(prctile(LCT,97.5))],'k','Facealpha',0.3,'LineStyle','none'); hold on
% plot([0:35],MLECT,'k','LineWidth',2);
% 
% patch([0:35 35:-1:0],[prctile(LNSCT,2.5) flip(prctile(LNSCT,97.5))],'r','Facealpha',0.3,'LineStyle','none'); hold on
% plot([0:35],MLECTNS,'r','LineWidth',2);
% 
% xlim([0 35]);
% xlabel('Days from exposure to isolation','Fontsize',18)
% yh=ylabel('Expected probability','Fontsize',18);
% box off;
% set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[0:6:42],'Xminortick','on');
% text(yh.Extent(1),max(ylim),'C','Fontsize',32,'FontWeight','bold');
% 
% subplot('Position',[0.397666849348435,0.17,0.588677688466691,0.341162790697676]); 
% for ii=1:4
%    patch(ii+[-0.35 -0.35 0.35 0.35],[prctile(LCT(:,[3.*ii+1]),[2.5 97.5 97.5 2.5])].*(Inc(end,2)+IncO(end,2)),'k','Facealpha',0.3,'LineStyle','none'); hold on
%    plot(ii+linspace(-0.35,0.35,101),MLECT(3.*ii+1).*(Inc(end,2)+IncO(end,2)).*ones(101,1),'k','LineWidth',2); hold on   
%    
%    patch(ii+[-0.35 -0.35 0.35 0.35],[prctile(LNSCT(:,[3.*ii+1]),[2.5 97.5 97.5 2.5])].*(Inc(end,2)+IncO(end,2)),'r','Facealpha',0.3,'LineStyle','none'); hold on
%    plot(ii+linspace(-0.35,0.35,101),MLECTNS(3.*ii+1).*(Inc(end,2)+IncO(end,2)).*ones(101,1),'r','LineWidth',2); hold on
% end
% xh=xlabel('Days from exposure to isolation','Fontsize',18);
% xh.Position=[2.500001907348632,-37.18,-1];
% yh=ylabel({'Expected number of','cases missed'},'Fontsize',18);
% box off;
% 
% ylim([0 300]);
% xlim([0.5 4.5]);
% set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1 2 3 4],'XTickLabel',{'3','6','9','12'},'YTick',[0:25:300]);
% xtickangle(45);
% text(yh.Extent(1),max(ylim),'D','Fontsize',32,'FontWeight','bold');