load('Incubation.mat','I')
test=fitdist(I','NegativeBinomial');
CV=test.ParameterCovariance;
G=test.ParameterValues;
save('Par_NBin.mat','G','CV');

clear;
load('Serial.mat','S')
test=fitdist(S','NegativeBinomial');
CVM=test.ParameterCovariance;
M=test.ParameterValues;
save('Par_NBin_Serial.mat','M','CVM');

clear;
