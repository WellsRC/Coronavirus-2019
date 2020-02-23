function [TH]=TimeHospitalDec31(S)
%% Distribution with onset prior up to December 31
dd=[13.9965508775716,1.76371520889767];
TH=round(wblrnd(dd(1),dd(2),S,1));

end
