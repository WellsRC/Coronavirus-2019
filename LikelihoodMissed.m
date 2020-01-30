function L = LikelihoodMissed(pt)
%LIKELIHOODMISSED Calculates the probability that the individual will be
%missed in airport screening

%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%
% pt - probability of travel
% s - duration of incubation period

%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%
% L - likelihood of travelling during the days before they become
% symptomatic

%%%%%%%%%%%%%%%%%%%%%%
% Calculation
%%%%%%%%%%%%%%%%%%%%%%
L=1-prod(1-pt); 
end

