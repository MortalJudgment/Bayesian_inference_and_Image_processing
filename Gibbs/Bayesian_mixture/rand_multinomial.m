function [distn, rndexps] = rand_multinomial(n,probs)
%----------------------------------------------
% [distn, rnsmult] = rand_multinomial(n, probs)
% Multinomial parameters n and probs (of m cells not necessarily summing to 1)
% distn - realization of rand_multinomial
% rndexps - random experiments, n of those with 1 in one of m cells
probs = probs(:)/sum(probs); % normalize
m = size(probs,1);
rndexps = zeros(m,n); %initialize random experiments
for j = 1:n  % assign 1 to one of m cells
    cumsumm=0;
    urand=rand;
  for i = 1:m
  cumsumm = cumsumm + probs(i,1);
    if  urand < cumsumm
        rndexps(i,j)=1;
        break
    end
  end
end
if n==1
    distn = rndexps'; %if n=1 take rndmult as realization
else
    distn = sum(rndexps'); %if n>1 sum over the experiments
end
% Bayes Stat at GaTech, June 2004------------------------------------------- 

