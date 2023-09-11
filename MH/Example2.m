close all; clear all;

% Inline functions
strg='1./(1+x.^2)';
cauchy=inline(strg,'x');

strg='1/sig*exp(-0.5*((x-mu)/sig).^2)';
norm=inline(strg,'x','mu','sig');

    
n=10000; %number of iterations
burn=500; %burn-in
sig=1;                  
x=zeros(1,n); %allocate memory for x      

%generate starting point from a standard normal distribution (To change)    
x(1)=12;
%x(1)=0.5;

for i=2:n  
    y=randn(1);% candidate
    u=rand(1);
    alpha=min([1, cauchy(y)*norm(x(i-1),0,sig)/(cauchy(x(i-1))*norm(y,0,sig))]); 
    if u<alpha; x(i)=y; 
    else x(i)=x(i-1); 
    end
end

figure;
plot(x(burn+1:n)) 
xlabel('iterations');ylabel('x');title('Post burn-in trajectory')

figure;
[N, xx]=hist(x(burn+1:n),[-20:0.5:20]); 
bar(xx, N/((n-burn)*0.5)); hold on
plot(xx, 1./(pi.*(1+xx.^2)),'-r');
xlabel('x');title('MH algo: target=C(0,1), proposal=N(0,1), burn-in=500');
legend('Histog.','target=Cauchy');hold off
