close all
clear all

%%
%% Visualization of data
load coal
figure(1); 
    plot(year,y,'.' ); 
    xlabel('year');ylabel('Number of Disasters');
    title('Number of coal mine Disaters per year')
                               
%% ALGO

n = length(y); % nb of data
m = 1100;	  % total nb of iterations

% Hyperparameters that are fixed
a1 = 0.5;
a2 = 0.5;
c1 = 0;
c2 = 0;
d1 = 1;
d2 = 1;

theta = zeros(1,m);  %Allocation of memory
lambda = zeros(1,m);  %Allocation of memory
k = zeros(1,m);       %Allocation of memory
like = zeros(1,n); % Holds probabilities for k.

% Get starting points.
k(1) = unidrnd(n,1,1);  % Note that k will indicate an index to the year 
                        ...that corresponds to a hypothesized change-point.
theta(1) = 1;
lambda(1) = 1;
b1 = 1; 
b2 = 1;

% Start the Gibbs Sampler.
for i = 2:m 
   kk = k(i-1);
   
   % Get parameters for generating theta
   t = a1 + sum(y(1:kk)); 
   lam = kk + b1;  
   % Generate the variate for theta.
   theta(i) = gamrnd(t,1/lam,1,1); 
   
   % Get parameters for generating lambda.
   t = a2 + sum(y) - sum(y(1:kk));
   lam = n-kk+b2;
   % Generate the variate for lambda.
   lambda(i) = gamrnd(t,1/lam,1,1);
   
   % Generate the parameters b1 and b2.
   b1 = gamrnd(a1+c1,1/(theta(i)+d1),1,1);
   b2 = gamrnd(a2+c2,1/(lambda(i)+d2),1,1);
   
   % Now get the probabilities for k.
   for j = 1:n
      like(j) = exp((lambda(i)-theta(i))*j)*(theta(i)/lambda(i))^sum(y(1:j));
   end
   like = like/sum(like);
   k(i) = discreteinvrnd(like,1,1); 
end  


%%  VISUALIZATION

%%%%%%%%%%% Histogram of k's %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
    hist(k,112);
    title('Histogram of k')

%%%%%%%%%%%%%%%%%%% Histograms of theta and lambda %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
    subplot(1,2,1)
    hist(theta(101:1100)); 
    title('Histogram of \theta')
    subplot(1,2,2)
    hist(lambda(101:1100));
    title('Histogram of \lambda')
    
    
   
mode(k)
kk = k(101:1100);
l1 = length(kk);
l2 = length(find(kk==41)); % number of time k = 41
l2/l1   % probability P(k=41)

thetas = theta(101:1100);
lambdas = lambda(101:1100);
length(find(thetas - lambdas>0))

