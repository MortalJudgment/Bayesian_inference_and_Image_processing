close all;
clear all;
n = 200;        % sample size
y = 4 * randn(1,n) + 1;
 
%% 
NN = 10000;     % number of iterations
mus = [];
taus = [];
suma = sum(y);
mu = 0;         % set the parameters as prior means
tau = 2;        % 
for i = 1 : NN
  new_mu    = sqrt(1/(1+n*tau)) * randn + (tau * suma)/(1+n*tau);
  mu        = new_mu;
  par       = 1+1/2 * sum ( (y - mu).^2 );
  new_tau   = gamrnd(2 + n/2, 1/par, 1,1);
  tau       = new_tau;
  
  mus = [mus new_mu];
  taus = [taus new_tau];
end

%% PLOT
figure
burnin = 1000;
subplot(2,1,1)
hist(mus(burnin+1:NN), 100); title('Histograms of \mu (after the burn-in)')
subplot(2,1,2)
hist(taus(burnin+1:NN), 100); title('Histograms of \tau (after the burn-in)')

muhat = mean(mus(burnin+1:NN))
tauhat = mean(taus(burnin+1:NN))