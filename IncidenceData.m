function [IncF, IncOF]=IncidenceData

load('Incidence_Start_Dec292019.mat')
load('CumulativeIncidence.mat')
INDX=datenum('01-21-2020')-datenum('12-29-2019')+1;
TempI=[China(2:end)-China(1:end-1)]'; % Starts Jan. 21 as index starts the 20th
TempI=[[INDX:(INDX+length(TempI)-1)]' TempI];
for ii=1:length(Inc(:,1))
    ff=find(Inc2(:,1)==Inc(ii,1));
    if(~isempty(ff))
        Inc(ii,2)=max([Inc(ii,2) Inc2(ff,2)]);
    end
    ff=find(TempI(:,1)==Inc(ii,1));
    if(~isempty(ff))
        Inc(ii,2)=max([Inc(ii,2) TempI(ff,2)]);
    end
end
Inc=[Inc;TempI(end-3:end,:)];
TempI=[Other(2:end)-Other(1:end-1)]'; % Starts Jan. 21 as index starts the 20th
TempI=[[INDX:(INDX+length(TempI)-1)]' TempI];

for ii=1:length(IncO(:,1))
    ff=find(TempI(:,1)==IncO(ii,1));
    if(~isempty(ff))
        IncO(ii,2)=max([IncO(ii,2) TempI(ff,2)]);
    end
end
IncO=[IncO; TempI(end-1:end,:)];

IncF=zeros(30,2);
IncOF=zeros(30,2);
for ii=1:30
    ff=find(Inc(:,1)==ii);
    IncF(ii,1)=ii;
    if(~isempty(ff))
    	IncF(ii,2)=Inc(ff,2);
    end
    
    ff=find(IncO(:,1)==ii);
    IncOF(ii,1)=ii;
    if(~isempty(ff))
    	IncOF(ii,2)=IncO(ff,2);
    end
end
end