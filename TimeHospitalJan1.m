function [TH]=TimeHospitalJan1(S)
%% Distribution with onset prior up to December 31
dd=[10.1934554700002,2.61866762107576];
TH=round(wblrnd(dd(1),dd(2),S,1));

end
