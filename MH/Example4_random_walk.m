%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Target=N(0,1) and three proposal distributions:
    % a)q(.|x)=N(X,0.5)
    % b)q(.|x)=N(X,0.1)
    % c)q(.|x)=N(X,10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;

strg='1/sig*exp(-0.5*((x-mu)/sig).^2)';
norm=inline(strg,'x','mu','sig');

 
n=500;  % number of iterations for each chain.  
sig=[0.5, 0.1, 10]; % the variances
X=zeros(1,n,3); % Set up a 3D vector to store the samples.


X(1,1,:)=[ -10, 0, 0]; % get the starting values for the chains
                       
figure                                     
for run=1:3 
   for i=2:n
      y=X(1,i-1,run)+randn(1)*sig(run);
      u=rand(1);
      alpha=norm(y,0,1)/norm(X(1,i-1,run),0,1);
      if u<alpha; X(1,i,run)=y;
      else X(1,i,run)=X(1,i-1,run);
      end
   end
   
subplot(3,1,run)
hold on
plot(X(1,:,run));
ylabel('x');title(['chain ', num2str(run)])
plot([1,n],[2,2],'r',[1,n],[-2,-2],'r')     
end

