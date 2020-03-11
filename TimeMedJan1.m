function [TH]=TimeMedJan1(S)
%% Distribution with onset Jan 1 and after for first medical visit
dd=[4.96191211844973,1.21314450792566];
TH=round(wblrnd(dd(1),dd(2),S,1));
end
