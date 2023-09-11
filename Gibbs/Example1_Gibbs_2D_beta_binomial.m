%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;

%%
% Set up preliminaries.     
% Note: we use k for the chain length, because n is used for the number of trials in the binomial.
k = 1000;		% number of iterations
m = 500; 		% burn-in
a = 2;	        % a and b are parameters of the beta
b = 4;
n = 16;         % binomial parameter
x = zeros(1,k); 
y = zeros(1,k); 
%% 
% Pick a starting point.
x(1) = binornd(n,0.5,1,1); 
y(1) = betarnd(x(1)+a,n-x(1)+b,1,1); 
for i = 2:k
		x(i) = binornd(n,y(i-1),1,1);
		y(i) = betarnd(x(i)+a,n-x(i)+b,1,1); 
end

%% Estimated marginal
fhat = zeros(1,17);
for i = 1:17
		fhat(i) = mean(binopdf(i-1,n,y(500:k))); 
end
%%
%% true marginal
xx= 0:16; 
f = factorial(n)./factorial(n-xx)./ factorial(xx).*gamma(a+b).*gamma(xx+a)...
    .*gamma(n-xx+b)./gamma(a)./gamma(b)./gamma(a+b+n);
%%

%plot the two distributions 
subplot(1,2,1)
bar(1:17,fhat,1)
title('estimated marginal');
colormap([0.5,0.8,0.3])

subplot(1,2,2)
bar(1:17,f,1);
title('true marginal')
colormap([0.5,0.8,0.3])
%%