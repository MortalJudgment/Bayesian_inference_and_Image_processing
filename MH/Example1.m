

close all;
clear all;

%%%%%%%%%%%%  Inline functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Cauchy density (without constant terms).
strg='1./(1+x.^2)';
cauchy=inline(strg,'x');

% Gaussian density (without constant terms)
strg='1/sig*exp(-0.5*((x-mu)/sig).^2)';
norm=inline(strg,'x','mu','sig');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                  
n=10000 ; % number of iterations
burn=500; % burn-in
sig=2;                  
x=zeros(1,n);    %allocate memory for x             
     
x(1)=randn(1);  % starting point
for i=2:n     
   y=x(i-1)+sig*randn(1); % candidate
   u=rand(1); 
  alpha=min([1,cauchy(y)*norm(x(i-1),y,sig)/(cauchy(x(i-1))*norm(y,x(i-1),sig))]);
   if u<alpha; x(i)=y; % proposal is accepted  so x(t)=y
   else x(i)=x(i-1);   % proposal is refused,  and x(t)=x(t-1)
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%  PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(x(burn+1:n)) 
xlabel('iterations');ylabel('x');title('Post burn-in trajectory')

%[N,xx]=hist(Y,X) 
figure;
[N, xx]=hist(x(burn+1:n),[-20:0.5:20]); 
bar(xx, N/((n-burn)*0.5)); hold on
plot(xx, 1./(pi.*(1+xx.^2)),'-r');
xlabel('x');title('MH algo: target=C(0,1), proposition="adjusted"Normal, burn-in=500');
legend('Histog.','target=Cauchy')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
