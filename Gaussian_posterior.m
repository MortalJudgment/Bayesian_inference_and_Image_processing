clc; close all; clear all;
n = 1000;

% values of the hyperparameters of the prior on mu
mu_0 = 1; sigma_0 = 0.5;

% parameter value for the data to be simulated
mu = 0; sigma = 2;  % to sample X1~N(0,2^2)

% sampling the data
x = mu + sigma.*randn(1,n);

% Parameters for the posterior of mu
xbar = sum(x)/n;
lambda = (n*sigma_0^2)/(n*sigma_0^2 + sigma^2);
mu_n = lambda*xbar + (1-lambda)*mu_0;

sigma_n = (sigma^2*sigma_0^2)/(n*sigma_0^2+sigma^2);

% Subdivision of the x-axis
naxe = 1000;
axe = linspace(-3,3,naxe);

% Plotting
figure
plot(axe,normpdf(axe,mu_0,sigma_0^2),'r');
hold on
plot(axe,normpdf(axe,mu_n,sigma_n),'b');
title('Prior and posterior of \mu ')
legend('prior of \mu','posterior of \mu')

% Bayesian estimate
mu_n