clc;
clear;
%% Probabaility of Travel
load('Probability_Travel.mat')
w=exp(F)./sum(exp(F));
wc=cumsum(w);
r=rand(10^4,1);
spc=zeros(10^4,1);
for ii=1:10^4
f=find(r(ii)<=wc);
f=f(1);
spc(ii)=pc(f);
end
fprintf('Probability of travel: %5.4f (95%% CI %5.4f - %5.4f0) \n',[pc(F==max(F)) prctile(spc,[2.5 97.5])]);

%% Probability travel periods
load('Missed_Screening_TimetoIsolation.mat')
fprintf('Probability of travel during incubation: %2.1f (95%% CI %2.1f - %2.1f) \n',100.*[MLE prctile(L,[2.5 97.5])]);
fprintf('Probability of travel during incubation and time to medical vist (before Jan 1): %2.1f (95%% CI %2.1f - %2.1f) \n',100.*[MLENS prctile(LNS,[2.5 97.5])]);
fprintf('Probability of travel during incubation and time to medical vist (after Jan 1): %2.1f (95%% CI %2.1f - %2.1f) \n',100.*[MLENSA prctile(LNSA,[2.5 97.5])]);
fprintf('Probability of travel during incubation and time to hospital vist (before Jan 1): %2.1f (95%% CI %2.1f - %2.1f) \n',100.*[MLENSH prctile(LNSH,[2.5 97.5])]);
fprintf('Probability of travel during incubation and time to hospital vist (after Jan 1): %2.1f (95%% CI %2.1f - %2.1f) \n',100.*[MLENSHA prctile(LNSHA,[2.5 97.5])]);

%% Temporal
load('Daily_Prob_Expect.mat');