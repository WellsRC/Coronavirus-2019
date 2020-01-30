function L = LikelihoodVisitorMissed(dv,A,G)
%LIKELIHOODMISSED Calculates the probability that the individual will be
%missed in airport screening

%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%
% dv - duration of time being spent in epicenter
% A - the attack rate for the time being spent in th epi center
% G - the parameters for the gamma distribtuion
%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%
% L - likelihood of being asympotamtic a tthe time of screening

%%%%%%%%%%%%%%%%%%%%%%
% Calculation
%%%%%%%%%%%%%%%%%%%%%%
L=0;
for ii=1:dv % calcuate the probability for each day and sum
    if(ii==1)
        L=L+A(ii).*(1-nbinpdf(dv-ii,G(1),G(2)));
    else
        L=L+A(ii).*prod(1-A(1:(ii-1))).*(1-nbincdf(dv-ii,G(1),G(2)));
    end
end
end

