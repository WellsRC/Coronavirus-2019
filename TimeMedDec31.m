function [TH]=TimeMedDec31(S)
%% Distribution with onset prior up to December 31
dd=[5.90974219009578,1.01776338647239];
TH=round(wblrnd(dd(1),dd(2),S,1));

end
