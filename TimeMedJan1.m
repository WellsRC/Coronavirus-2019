function [TH]=TimeMedJan1(S)
%% Distribution with onset prior up to December 31
dd=[4.96191211844973,1.21314450792566];
TH=round(wblrnd(dd(1),dd(2),S,1));
end
