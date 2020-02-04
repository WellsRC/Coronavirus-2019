load('Weighted_Travel_Inubation');
load('Weight_Flights');
minE=-22;
maxE=57;
EsT=zeros(10,2);
for ii=1:length(FC)
    tf = strcmp({FC{ii,1}},{FlightAll{:,1}});
    wt=FlightAll{tf,2};
    PP=zeros(size(MLE(1,:)))';
    for jj=1:length(PP)
       PP(jj)=pchip(w,MLE(:,jj),wt); 
    end
    EsT(ii,1)=[minE:maxE]*PP./sum(PP);
    EsT(ii,2)=FC{ii,2};
end
f1=fopen('DatesofExportation.txt','w');
startDateofSim = datenum('12-06-2019');% Start date
for ii=1:length(FC)
    fprintf(f1,[ FC{ii,1} ',' datestr([startDateofSim+(EsT(ii,1)-1)],'mmm-dd-yy') ',' datestr([startDateofSim+(EsT(ii,2)-1)],'mmm-dd-yy') '\n']);
    fprintf([ FC{ii,1} ',' datestr([startDateofSim+(EsT(ii,1)-1)],'mmm-dd-yy') ',' datestr([startDateofSim+(EsT(ii,2)-1)],'mmm-dd-yy') '\n'])
end

fclose(f1);
