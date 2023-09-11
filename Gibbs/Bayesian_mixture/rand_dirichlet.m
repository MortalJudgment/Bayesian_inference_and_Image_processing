function drand = rand_dirichlet(a,n)
% function drand = rand_dirichlet(a)
% a - vector of parameters.
%---------------------------------------------------
a=a(:);
m=size(a,1);
a1=zeros(m,n);
for i = 1:m
    a1(i,:)=rand_gamma(a(i,1),1,1,n);
end
for i=1:m
drand(i, 1:n )= a1(i, 1:n ) ./ sum(a1);
end

    

