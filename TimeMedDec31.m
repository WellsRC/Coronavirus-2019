function [TH]=TimeMedDec31(S)
%% Distribution with onset prior to and including December 31 for first medical visit
dd=[5.90974219009578,1.01776338647239];
TH=round(wblrnd(dd(1),dd(2),S,1));

end
