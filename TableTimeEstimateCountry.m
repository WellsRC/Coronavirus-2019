clc;
load('Weight_Flights');
minE=-22;
maxE=57;
EsT=zeros(10,3);
    load(['Weighted_Travel_Infectious_Country.mat']);
    for ii=1:length(FC)
        tf = strcmp({FC{ii,1}},{FlightAll{:,1}});
        wt=FlightAll{tf,2};
        f=find(w==wt);
        PP=MLE(f,:)';
        PP=PP(2:end)-PP(1:(end-1));
        PP=PP./sum(PP);
        EsT(ii,1)=([(minE+1):maxE]*PP);
        EsT(ii,2)=sqrt(([(minE+1):maxE]).^2*PP-(EsT(ii,1)).^2);
        EsT(ii,3)=FC{ii,2};
    end
f1=fopen('DatesofExportation.txt','w');
startDateofSim = datenum('12-06-2019');% Start date
for ii=1:length(FC)
    fprintf(f1,[ FC{ii,1} ',' datestr([startDateofSim+(EsT(ii,1)-1)],'yyyy-mm-dd') ', %2.1f,' datestr([startDateofSim+(EsT(ii,3)-1)],'yyyy-mm-dd') ' \n'],EsT(ii,2));
    fprintf([ FC{ii,1} ',' datestr([startDateofSim+(EsT(ii,1)-1)],'yyyy-mm-dd') ', %2.1f,' datestr([startDateofSim+(EsT(ii,3)-1)],'yyyy-mm-dd') ' \n'],EsT(ii,2));
end

fclose(f1);
