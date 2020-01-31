function [TH]=TimeHospitalDec31(S)
%% Distribution with onset prior up to December 31
load('TimetoHospitaluptoDec31','C');
r=rand(S,1);
TH=zeros(S,1);
for ii=1:S
    f=find(r(ii)<C);
    f=f(1);
    TH(ii)=f-1; % Subtract since the first position is zero
end

end
