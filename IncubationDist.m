function [I,pdf] = IncubationDist(m,S)
 % 95% CI of m is 4.1 to 7.0
    v=15.289588142870238; % produces the 95th percentile of 12.5 for a mean of 5.2 i.e logninv(0.95,mu,sigma)=12.5;
    mu = log((m^2)/sqrt(v+m^2));
    sigma = sqrt(log(v/(m^2)+1));
    I=[];
    if(S>0)
        I=round(lognrnd(mu,sigma,S,1));
    end
    pdf=[0:42];
    pdf(1)=logncdf(0.5,mu,sigma);
    pdf(2:43)=logncdf(0.5+[1:42],mu,sigma)-logncdf(0.5+[0:41],mu,sigma);
    
end

