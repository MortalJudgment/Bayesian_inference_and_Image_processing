clc
clear all
close all

% sampling the data
lambda1_true = 20;
lambda2_true = 13;
l1 = 100;
l2 = 100;

y = [poissrnd(lambda1_true,1,l1) poissrnd(lambda2_true,1,l2)];

%
Nmc = 5000; % number of  iterations
alpha = 2;

% initialization
lambda1 = 1;
lambda2 = 1;
beta = 1;


%%
for m_compt = 1:Nmc
    
    lambda1 = gen_gamma(sum(y(1:l1))+alpha, beta+l1);
    Tlambda1(m_compt) = lambda1;
    
    lambda2 = gen_gamma(sum(y((l1+1):end))+alpha, beta+l2);
    Tlambda2(m_compt) = lambda2;
    
    beta = gen_gamma(2*alpha, lambda1+lambda2);
    Tbeta(m_compt) = beta;
end

%%  Plotting
%
Nbi = 250; % burn-in
Nr = Nmc - Nbi;
Tlambda1 = Tlambda1(Nbi+1:end);
Tlambda2 = Tlambda2(Nbi+1:end);
Tbeta = Tbeta(Nbi+1:end);

beta_MMSE = mean(Tbeta);


figure(1)
    % estimated posterior
    pas1 = 0.1;
    min1 = 16;
    max1 = 24;
    vect1 = min1:pas1:max1;
    h1 = hist(Tlambda1,vect1)/(Nr*pas1);
    plot(vect1,h1)
    % true posterior
    pdf1 = gampdf(vect1,sum(y(1:l1))+alpha,1/(beta_MMSE+l1));
    hold on
    plot(vect1,pdf1,'r')
    xlabel('\lambda_1')
    ylabel('f(\lambda_1|\alpha, \beta,x)')
    legend('Estimated',' True')
    hold off
    axis([min1 max1 0 1.1*max(pdf1)])
    title('FCD of \lambda_1 (\lambda_1^*=20)')


figure(2)
% estimated posterior
    pas2 = 0.1;
    min2 = 10;
    max2 = 16;
    vect2 = min2:pas2:max2;
    h2 = hist(Tlambda2,vect2)/(Nr*pas2);
    plot(vect2,h2)
    % true posterior
    pdf2 = gampdf(vect2,sum(y(l1+1:end))+alpha,1/(beta_MMSE+l2));
    hold on
    plot(vect2,pdf2,'r')
    hold off
    xlabel('\lambda_2')
    ylabel('f(\lambda_2|\alpha,\beta,x)');
    legend('Simulated','True')
    axis([min2 max2 0 1.1*max(pdf2)])
    title('FCD  of \lambda_2 (\lambda_2^*=13)')


figure(3)
    % estimated posterior
    pas3 = 0.01;
    min3 = 0;
    max3 = 0.5;
    vect3 = min3:pas3:max3;
    h3 = hist(Tbeta,vect3)/(Nr*pas3);
    plot(vect3,h3)
    % true posterior
    pdf3 = gampdf(vect3,2*alpha,1/(lambda1_true+lambda2_true));
    hold on
    plot(vect3,pdf3,'r')
    hold off
    xlabel('\beta')
    ylabel('f(\beta|\lambda_1,\lambda_2,Y)');
    legend('Simulated','True')
    axis([min3 max3 0 1.1*max(pdf3)])
    title('FCD of \beta')
%%
