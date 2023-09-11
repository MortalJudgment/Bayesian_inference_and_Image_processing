close all
clear all

y=zeros(1,13);
lambda=1;
for i = 1 : 13
    while(y(i)<4 )
        y(i) = poissrnd(lambda,1,1);
    end
end
niter=2000;
lambdas=zeros(1, niter); % for the  non-Rao-Blackwell estimates
lambdasest=zeros(1, niter); % for the  Rao-Blackwell estimates

for j = 1 :niter
   lambda= gamrnd(313+sum(y), 1/360, 1,1);
   y=zeros(1,13);
   for i = 1 : 13
      while(y(i)<4 )
        y(i) = poissrnd(lambda,1,1);
       end
   end
   lambdahat=(313+sum(y))/360;
   lambdasest(j)=lambdahat;
   lambdas(j)=lambda;
end


% Convergence
lambdabarRB=zeros(1,niter);
lambdabar=zeros(1,niter);
for i = 1:niter 
		lambdabarRB(i) = mean(lambdasest(1:i));
        lambdabar(i) = mean(lambdas(1:i));
end
figure(1);plot(1:niter,lambdabar);xlabel('iterations');ylabel('average'); 
hold on
plot(1:niter,lambdabarRB, 'r')
legend('usual','R-B')
title('Evolution of the R-B empirical average')


% Estimates and histograms
lambdaestRB=mean(lambdasest(1000:2000))  %Rao-Blackwell estimate
lambdaest=mean(lambdas(1000:2000)) % Usual estimate

figure; hist(lambdasest(1000:2000),6)
figure; hist(lambdas(1000:2000),6)