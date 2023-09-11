close all;
clear all;

% Inline function

strg='1./(1+x.^2)';
cauchy=inline(strg,'x');

              
n=10000;%number of iterations
burn=500;%burn-in               
x=zeros(1,n); %allocate memory for x      

% starting point      
x(1)=12;
for i=2:n  
    y=trnd(0.5);% candidate
    u=rand(1);  
    alpha=min([1, cauchy(y)*tpdf(x(i-1),0.5)/(cauchy(x(i-1))*tpdf(y,0.5))]);  
    if u<alpha; x(i)=y; 
    else x(i)=x(i-1); 
    end
end


figure;
plot(x(burn+1:n)) 
xlabel('iterations');ylabel('x');title('post burn-in trajectory')


figure;
[N, xx]=hist(x(burn+1:n),[-20:0.5:20]); 
bar(xx, N/((n-burn)*0.5)); hold on
plot(xx, 1./(pi.*(1+xx.^2)),'-r');
xlabel('x');title('MH Algo: target=C(0,1), proposition=T(0.5), burn-in=500');
legend('Histog.','target=Cauchy')
hold off