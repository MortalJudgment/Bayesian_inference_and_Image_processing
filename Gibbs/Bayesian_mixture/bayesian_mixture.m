%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
randn('seed',1) % sample standard normal distribution

% simulating the data
n = 400; % number of observations
us = rand(1,n); % sample uniform distribution
for i=1:n
    if us(i) < 0.1
        x(i) = 2 * randn +  1; % 0.1 from N(1,2^2)
        % If Y ~ N(0,1) then sigma*Y+mu ~ N(mu,sigma^2)
    elseif us(i) > 0.8
        x(i) = 2 * randn + 20; % 0.2 from N(20,2^2)
    else
        x(i) = 2 * randn +  9; % 0.7 from N(9,2^2)
    end
end
X = sort(x);

K=3; % number of components

% hyperparameters:
alpha  = 5*ones(1,K);
tau    = 12;
tau2   = tau^2;
gamma  = 10;
delta  = 30;
burn   = 1000;
N      = 6000;  % total number of iterations

% initialisation of parametres
sig2 = 5;
omega= 1./K*ones(1,K);
mu   = 10*ones(1,K);

s = ones(1,n);
classe = 1;
for i=1:n
    s(i)=classe;
    if(mod(i,ceil(n/K))==0), classe = classe+1;
    end
end

mus = [];
sigs = [];
ws = [];

h = waitbar(0, 'Simulation in progress');

for iter=1:N
    for i=1:n
        for j=1:K
            sij(i,j)= (s(i)==j);
        end
    end
    
    
    nj = sum(sij);
    a = alpha+nj;
    w = rand_dirichlet(a,1);
    
    for j=1:K
        sj(j) = sum(X(s==j));
    end
    
    var = 1./( nj/sig2 + 1/tau2);
    moy = (tau2.*sj)./(nj*tau2+sig2);
    mu  = moy+sqrt(var).*randn([1,K]);
    
    for j=1:K
        sm2j(j)=sum((X(s==j)-mu(j)).^2);
    end
    
    sig2 = 1./rand_gamma(n/2+gamma,delta+(1/2)*sum(sm2j),1,1);
    
    for i=1:n
        for j=1:K
            pr(i,j)=w(j)*exp(-1/(2*sig2)*(X(i)-mu(j))^2);
        end
    end
    
    tot=[];
    for j=1:K
        tot=[tot,sum(pr')'];
    end
    pr=pr./tot;
    
    for j=1:n
        [aa,bb]=rand_multinomial(1,pr(j,:));
        s(j)=find(bb==1);
    end
    
    mus = [mus;mu];
    sigs = [sigs,sig2];
    ws = [ws,w];
    waitbar(iter/N)
end


% Plotting



% Estimates
% Mean
muhat=[];
for k=1:K
    muhat(k) = mean(mus(burn+1:N,k));
end

%Variance
sighat = mean(sigs(burn+1:N));

%Weights
what = [];
for k=1:K
    what(k)=mean(ws(k,burn+1:N));
end

% Estimated and true densities
figure
histo(X,30,0,1);
mw=what
mmu=muhat
msig2=sighat

hold on
cee=[-5:0.1:25];
% Estimated
est=mw(1).*1./sqrt(2 * pi * msig2).*exp(-1/(2*msig2) * (cee - mmu(1)).^2)+...
    mw(2).*1./sqrt(2 * pi * msig2).*exp(-1/(2*msig2) * (cee - mmu(2)).^2)+...
    mw(3).*1./sqrt(2 * pi * msig2).*exp(-1/(2*msig2) * (cee - mmu(3)).^2);
plot2=plot(cee, est,'g-');

% True
theo = 0.1 * 1./sqrt(2 * pi * 4).*exp(-1/(2*4) * (cee - 1).^2)+...
    0.2 * 1./sqrt(2 * pi * 4).*exp(-1/(2*4) * (cee - 20).^2)+...
    0.7 * 1./sqrt(2 * pi * 4).*exp(-1/(2*4) * (cee - 9).^2);
plot3=plot(cee, theo,'r--');
legend([plot2 plot3], 'estimated','true','Location', 'NorthWest');
hold off
