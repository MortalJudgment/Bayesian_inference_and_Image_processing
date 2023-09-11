% The data
Failures = [5, 1, 5, 14, 3, 19, 1, 1, 4, 22];
Times = [94.32, 15.72, 62.88, 125.76, 5.24, 31.44, 1.05, 1.05, 2.10, 10.48];

Nx = length(Failures);

%Hyperparameters
alpha = 1.8; gamma = 0.01; delta = 1;

Nsim = 5000; % Number of iterations

lambda = zeros(Nsim, Nx);
beta = zeros(1, Nsim);

%Initialization
beta(1) = gamrnd(gamma, 1/delta, 1,1);

%Gibbs
for j = 2 : Nsim
    for i = 1 : Nx
        par1 = Times(i) + beta(j-1);
        lambda(j,i) = gamrnd(Failures(i)+alpha, 1/par1, 1,1);
        par2 = delta + sum(lambda(j,:));
        beta(j) = gamrnd(gamma+Nx*alpha, 1/par2, 1,1);
    end
end


%Plotting and estimating
burnin = 500;

figure(1)
hist(beta(burnin+1:Nsim)); 
title('Histogram of \beta')

lambda1 = lambda(:,1);
figure(2)
hist(lambda1(burnin+1:Nsim)); 
title('Histogram  of \lambda_1')

lambda2 = lambda(:,2);
figure(3)
hist(lambda2(burnin+1:Nsim)); 
title('Histogram of \lambda_2')

