function [R,P] = SampleDistShape(G,CV,N)

P=randn(1.5.*N,1).*sqrt(CV(2,2))+G(2);
t1=(randn(1.5.*N,1)+CV(1,2).*P./sqrt(CV(2,2))-CV(1,2).*G(2)./sqrt(CV(2,2)));
sc=std(t1);
R=(sqrt(CV(1,1))./sc).*t1+G(1);

R=R(P>=0);
P=P(P>=0);
R=R(P<=1);
P=P(P<=1);

P=P(R>0);
R=R(R>0);

P=P(1:N);
R=R(1:N);
end

