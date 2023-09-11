% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
%%
n       = 5000;  % nb of iterations
m       = 500;   % burn in
rho     = 0.5; 
xgibbs  = zeros(2,n);   % xgibbs is a 2*n matrix, the first row of the matrix
                        % will contain the generated x's while the second will contain the y's
y       = [0;0];
sig     = sqrt(1-rho^2);%variance

%% ALGO

% Initial point.
xgibbs(:,1) = [10,10]; 

% Start the chain.
for i=2:n
     %Iteration i: sampling of x1(iteration i) given x2(iteration i-1) 
     mu = y(1)+rho*(xgibbs(2,i-1)-y(2));    % mean 
     xgibbs(1,i) = mu+sig*randn(1);         % sampling
     
     %Iteration i: sampling of x2(iteration i) given x1(iteration i)
     mu = y(2)+rho*(xgibbs(1,i)-y(1));
     xgibbs(2,i) = mu+sig*randn(1);  
end


%% Vizualisation
figure(1); 
plot(xgibbs(1,m+1:n), xgibbs(2,m+1:n), '.' );
xlabel('x_1');ylabel('x_2');
title(['Gibbs for a bivariate Normal, burn-in=' num2str(m) ', nb iterations=' num2str(n) ])  

figure (2);
subplot(2,1,1); hist(xgibbs(1,m+1:n),100);
title(['Histogram of x_1, burn-in=' num2str(m) ', nb iterations=' num2str(n) ])
subplot(2,1,2); hist(xgibbs(2,m+1:n),100);
title(['Histogram of  x_2, burn-in=' num2str(m) ', nb iterations=' num2str(n) ]);

x1 = xgibbs(1,:);
x2 = xgibbs(2,:);
xbar_1 = zeros(1,n);
xbar_2 = zeros(1,n);

for i = 1:n
    xbar_1(i) = mean(x1(1:i));
    xbar_2(i) = mean(x2(1:i));
end
figure;
plot(1:1000, xbar_1(1:1000))
hold on 
plot(1:1000, xbar_2(1:1000))
legend('mean of x','mean of y')

muhat = mean(xgibbs(:,m+1:n),2)
sigmahat = cov(xgibbs(:,m+1:n)')