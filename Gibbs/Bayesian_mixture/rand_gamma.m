function r = rand_gamma(a,b,m,n)
%   r = randgamm(a,b,m,n) returns a matrix of size (m,n)
%   with entries distributed as Gamma(a,b).
%     a == shape parameter
%     b == 1/scale parameter
%   ----------------------
%     mean = a /b; var = a /b^2;
%   ----------------------
%   density f(x) = x^(a-1)*(b^a) / Gamma(a) * exp(- b* x);
%   For example, Chi2_df is Gamma(df/2, 1/2)
if nargin < 2, 
   error('Requires at least two input arguments.'); 
end


if nargin == 2
   [errorcode rows columns] = rndcheck(2,2,a,b);
end

if nargin == 3
   [errorcode rows columns] = rndcheck(3,2,a,b,m);
end

if nargin == 4
   [errorcode rows columns] = rndcheck(4,2,a,b,m,n);
end

if errorcode > 0
   error('Size information is inconsistent.');
end

% Initialize R to zero.
lth = rows*columns;
r = zeros(lth,1);
a = a(:); b = 1./b(:);

scalara = (length(a) == 1);
if scalara 
   a = a*ones(lth,1);
end

scalarb = (length(b) == 1);
if scalarb 
   b = b*ones(lth,1);
end

% If a == 1, then gamma is exponential. (Devroye, page 405).
k = find(a == 1);
if any(k)
   r(k) = -b(k) .* log(rand(size(k)));
end 


k = find(a < 1 & a > 0);
% (Devroye, page 418 Johnk's generator)
if any(k)
	c = zeros(lth,1);
	d = zeros(lth,1);
   c(k) = 1 ./ a(k);
   d(k) = 1 ./ (1 - a(k));
   accept = k;
   while(length(accept)>0),
      u = rand(size(accept));
      v = rand(size(accept));
      x = u .^ c(accept);
      y = v .^ d(accept);
      k1 = find((x + y) <= 1); 
      if ~isempty(k1)
	      e = -log(rand(size(k1))); 
   	   r(accept(k1)) = e .* x(k1) ./ (x(k1) + y(k1));
         accept(k1) = [];
      end
   end
   r(k) = r(k) .* b(k);
end

% Use a rejection method for A > 1.
k = find(a > 1);
% (Devroye, page 410 Best's algorithm)
bb = zeros(size(a));
c  = bb;
if any(k)
   bb(k) = a(k) - 1;
   c(k) = 3 * a(k) - 3/4;
   accept = k; 
   count = 1;
   while(length(accept)>0)
      m = length(accept);
      u = rand(m,1);
      v = rand(m,1);
      w = u .* (1 - u);
      y = sqrt(c(accept) ./ w) .* (u - 0.5);
      x = bb(accept) + y;
      k1 = find(x >= 0);
      if ~isempty(k1)
	 z = 64 * (w .^ 3) .* (v .^ 2);
         k2 = (z(k1) <= (1 - 2 * (y(k1) .^2) ./ x(k1)));
         k3 = find(k2);
         k3 = k1(k3);
	 r(accept(k3)) = x(k3); 
         k4 = find(~k2);
         k4 = k1(k4);
	 k5 = k4(find(log(z(k4)) <= (2*(bb(accept(k4)).*log(x(k4)./bb(accept(k4)))-y(k4)))));
         r(accept(k5)) = x(k5);
         omit = [k3; k5];
         accept(omit) = [];
      end
   end
   r(k) = r(k) .* b(k);
end

% Return NaN if b is not positive.
if any(any(b <= 0));
   if prod(size(b) == 1)
      tmp = NaN;
      r = tmp(ones(rows,columns));
   else
      k = find(b <= 0);
      tmp = NaN;
      r(k) = tmp(ones(size(k)));
   end
end

% Return NaN if a is not positive.
if any(any(a <= 0))
   if prod(size(a) == 1)
      tmp = NaN;
      r = tmp(ones(rows,columns));
   else
      k = find(a <= 0);
      tmp = NaN;
      r(k) = tmp(ones(size(k)));
   end
end

r = reshape(r,rows,columns);
