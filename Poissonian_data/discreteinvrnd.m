function X = discreteinvrnd(p,m,n)
%generates (m*n) random numbers from any discrete distribution with proba
% mass vector given by p (based on an inversion method)
X=zeros(m,n);%preallocate memory for X
for i=1:m*n
    u=rand;
    I=find(u<cumsum(p));
    X(i)=min(I);
end
