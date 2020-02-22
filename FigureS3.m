load('Probability_Travel_Infection_6733_Mulit.mat')
close all;

startDateofSim = datenum('12-06-2019');% Start date
XTL=datestr([startDateofSim+[36:2:72]-1],'mm-dd-yy'); % Need to subtract one as index 1 is the start date Dec 6 2019


figure('units','normalized','outerposition',[0 0 1 1]);

subplot('Position',[0.055,0.165,0.283865546218487,0.3]); 

contourf([36:72],pc,F','LineStyle','none');
xlabel('Date data included up to','Fontsize',18);
yh=ylabel('Probability of travel per day','Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[36:2:72],'XTickLabel',XTL,'Xminortick','on','Yminortick','on');
xtickangle(45);
y=colorbar;


yhC=ylabel(y,{'log-likelihood'});
yhC.Rotation=270;
yhC.Position=[4.357460260391235,-200.7109392296664,0];
y.Position=[0.344432011129055,0.164212573539994,0.011204481792717,0.29971180782907];
text(yh.Extent(1),0.032,'A','Fontsize',32,'FontWeight','bold');
colormap('gray');

subplot('Position',[0.454,0.165,0.283865546218487,0.3]); 
MM=zeros(37,3);
r=rand(10000,1);
spc=zeros(10000,1);
for ii=1:37
    MM(ii,1)=find(F(ii,:)==max(F(ii,:)));
    w=exp(F(ii,:))./sum(exp(F(ii,:)));
    wc=cumsum(w);
    for jj=1:10000
        f=find(r(jj)<=wc);
        f=f(1);
        spc(jj)=pc(f);
    end
    MM(ii,2:3)=prctile(spc,[2.5 97.5]);
end
plot([36:72],pc(MM(:,1)),'k','LineWidth',2); hold on
patch([36:72 72:-1:36],[(MM(:,2))' flip((MM(:,3)'))],'k','Facealpha',0.35,'LineStyle','none');
xlabel('Date data included up to','Fontsize',18);
ylabel('Probability of travel per day','Fontsize',18);
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[36:2:72],'XTickLabel',XTL,'Xminortick','on','Yminortick','on');
xtickangle(45);
xlim([36 72]);
ylim([0 0.02]);

text(yh.Extent(1),0.0214,'B','Fontsize',32,'FontWeight','bold');