function L = LikelihoodScreen(pt,s)
%LIKELIHOODScreen Calculates the probability that the individual will be
%caught in airport screening

%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%
% pt - probability of travel
% s - duration of incubation period

%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%
% L - likelihood of being symptomatic on the day of travel

%%%%%%%%%%%%%%%%%%%%%%
% Calculation
%%%%%%%%%%%%%%%%%%%%%%
L=pt(1-pt).^s; 
end

