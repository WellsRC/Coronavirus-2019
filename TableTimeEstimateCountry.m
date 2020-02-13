% clc;

EsT=zeros(10,5);
    load(['Weighted_Travel_Infectious_Country.mat']);
    load('Weight_Flights');

minE=-22;
maxE=68;
    for ii=1:length(FC)
        PP=MLE(ii,:)';
        PP=PP(2:end)-PP(1:(end-1));
        PP=PP./sum(PP);
        tt=([(minE+1):maxE]);
        EsT(ii,4)=tt(PP==max(PP));
        EsT(ii,3)=([(minE+1):maxE]*PP);
        EsT(ii,5)=sqrt(([(minE+1):maxE]).^2*PP-(EsT(ii,3)).^2);
        EsT(ii,1)=FC{ii,2};
        EsT(ii,2)=wt;
    end
f1=fopen('DatesofExportation.txt','w');
startDateofSim = datenum('12-06-2019');% Start date
for ii=1:length(FC)
    fprintf(f1,[ FC{ii,1} ',%3.2E ,' datestr([startDateofSim+(EsT(ii,1)-1)],'yyyy-mm-dd') ', ' datestr([startDateofSim+(EsT(ii,4)-1)],'yyyy-mm-dd') ',' datestr([startDateofSim+(EsT(ii,3)-1)],'yyyy-mm-dd') ', %2.1f \n'],[EsT(ii,2) EsT(ii,5)]);
    fprintf([ FC{ii,1} ',%3.2E ,' datestr([startDateofSim+(EsT(ii,1)-1)],'yyyy-mm-dd') ' ,' datestr([startDateofSim+(EsT(ii,4)-1)],'yyyy-mm-dd') ',' datestr([startDateofSim+(EsT(ii,3)-1)],'yyyy-mm-dd') ', %2.1f \n'],[EsT(ii,2) EsT(ii,5)]);
end

fclose(f1);
