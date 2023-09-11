             % udeconv.m --- Unsupervised Wiener-Hunt deconvolution

% Copyright  (C)  2011  François Orieux <francois.orieux@gmail.com>

% Version: 1.1
% Keywords: Deconvolution, Unsupervised, MCMC, Wiener-Hunt
% Author: François Orieux <francois.orieux@gmail.com>
% Maintainer: François Orieux <francois.orieux@gmail.com>
% URL: http://www.lss.supelec.fr/perso/orieux_francois/index.html

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:

% Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.

% Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.

% Neither the name of François Orieux nor the names of its
% contributors may be used to endorse or promote products derived from
% this software without specific prior written permission.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
% COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% Code:

function [xEap, gnChain, gxChain, varargout] = udeconv(data, ir, varargin)
% udeconv - Unsupervised Wiener-Hunt deconvolution
%
% xEap = udeconv(data, ir, ...)
% [xEap, gnChain, gxChain] = udeconv(data, ir, ...)
%
% return the deconvolution of 'data' by 'ir'. The algorithm is a
% stochastic iterative process (Gibbs sampler) that allow automatic
% tuning of regularization parameter, see reference below. There is no
% specific constraints on the number of dimension.
%
% The call [xEap, gnChain, gxChain, xStd] = udeconv(...) allow to
% compute the diagonal of the covariance matrix around xEap with the
% cost of an fft at each iteration.
%
% If you use this work, add a citation of the reference below.
%
% Compatible with octave.
%
% PARAMETERS
%
% data -- the data
%
% ir -- the impulsionnal response
%
% OPTIONNALS
%
% Optionnals argument are in the form (..., 'key', val, ...).
%
% 'criterion', val -- if the difference between two successive estimate is
% less than this value, stop the algorithm. Default is 1e-4.
%
% 'burnin', val -- number of iteration to remove at the beginning of the
% chain to compute the mean of the image. Default is typicaly 30.
%
% 'maxIter', val -- maximum number of iteration. Default is 150.
%
% OUTPUTS
%
% xEap -- the estimated result
%
% gnChain, gxChain -- the MCMC chain of the regularisation
% parameters. See reference below.
%
% xStd -- is the standart deviation around the estimate
%
% FUNCTION CALL
%
% xEap = udeconv(data, ir)
%
% [xEap gnChain, gxChain] = udeconv(...)
%
% [xEap gnChain, gxChain, xStd] = udeconv(...)
%
% REFERENCE
%
% François Orieux, Jean-François Giovannelli, and Thomas Rodet,
% "Bayesian estimation of regularization and point spread function
% parameters for Wiener-Hunt deconvolution," J. Opt. Soc. Am. A 27,
% 1593-1607 (2010)
%
% http://www.opticsinfobase.org/josaa/abstract.cfm?URI=josaa-27-7-1593

  criterion = 1e-4;
  burnin = 30;
  maxIter = 150;

  for iargin = 3:2:nargin
    if strcmp(varargin{iargin-2}, 'criterion')
      criterion = varargin{iargin-1};
    elseif strcmp(varargin{iargin-2}, 'burnin')
      burnin = varargin{iargin-1};
    elseif strcmp(varargin{iargin-2}, 'maxIter')
      maxIter = varargin{iargin-1};
    else
      disp('Option not available')
    end
  end

  % The mean of the object
  circXeap = zeros(size(data));
  % The previous computed mean in the iterative loop
  previousCircXeap = zeros(size(data));

  xStd = zeros(size(data));
  previousXstd = zeros(size(data));

  tf = ir2tf(ir, size(data));
  hpFilter = ir2tf([0 -1 0; -1 4 -1; 0 -1 0], size(data));

  % Difference between two succesive mean
  delta = NaN;

  % Initial state of the chain
  gxSample = 1; gnSample = 1;
  nPixels = numel(data);

  % Parameter of the hyperparameter law. The following value correspond
  % to Jeffery's prior for the hyper parameter. See reference.
  alphaN = 0; betaNbar = 0; % = 1/betaB
  alphaX = 0; betaXbar = 0; % = 1/beta1

  % The correlation of the object in Fourier space (if size is big,
  % this can reduce computation time in the loop)
  ahpFilter2 = abs(hpFilter).^2;
  atf2 = abs(tf).^2;

  data = ufft(data);

  % Gibbs sampling
  for iteration = 1:maxIter
    %% Sample of Eq. 27 p(circX^k | gn^k-1, gx^k-1, y).

    % white complex gaussian noise (no respect of hermitian symmetry)
    circXcenteredSample = sqrt(0.5)*(randn(size(data)) + ...
                                     sqrt(-1)*randn(size(data)));

    % ponderation (correlation in direct space)
    inverseVariance = gnSample*atf2 + gxSample*ahpFilter2; % Eq. 29
    circXcenteredSample = circXcenteredSample./sqrt(inverseVariance);

    % mean Eq. 30 (RLS for fixed gn, gamma0 and gamma1 ...)
    gain = gnSample*conj(tf)./inverseVariance;
    circXmean = gain.*data;

    % Sample of X in Fourier space
    circXsample = circXcenteredSample + circXmean;

    %% Sample of Eq. 31 p(gn | x^k, gx^k, y)
    likelihood = sum(abs(data(:) - circXsample(:).*tf(:)).^2);
    gnSample = gamrnd(alphaN + nPixels/2, 1/(betaNbar + likelihood/2));
    gnChain(iteration) = gnSample;

    %% Sample of Eq. 31 p(gx | x^k, gn^k-1, y)
    smoothness = sum(abs(circXsample(:).*hpFilter(:)).^2);
    gxSample = gamrnd (alphaX + (nPixels - 1)/2, 1/(betaXbar + ...
                                                    smoothness/2));
    gxChain(iteration) = gxSample;

    %% Current empirical average
    if (iteration > burnin)
      circXeap = previousCircXeap + circXsample;

      if (iteration > (burnin + 1))
        norm = sum(abs(circXeap(:)))/(iteration - burnin);
        current = circXeap/(iteration - burnin);
        previous = previousCircXeap/(iteration - burnin - 1);

        delta = sum(abs(current(:) - previous(:)))/norm;
      end

      previousCircXeap = circXeap;

      if (nargout == 4)
        xStd = previousXstd + real(uifft(circXsample)).^2;
        previousXstd = xStd;
      end

    end

    %% Stop of the algorithm
    if (delta < criterion)
      circXeap = circXeap/(iteration - burnin);
      xEap = real(uifft(circXeap));

      if (nargout == 4)
        xStd = xStd/(iteration - burnin);
        xStd = sqrt(xStd - xEap.^2);
        varargout = {xStd};
      end

      return
    end
  end

  %% Empirical average \approx EAP Eq. 44
  circXeap = circXeap/(maxIter - burnin);
  xEap = real(uifft(circXeap));

  if (nargout == 4)
    xStd = xStd/(iteration - burnin);
    xStd = sqrt(xStd - xEap.^2);
    varargout = {xStd};
  end

end

function y = ufft(x, varargin)
% UFFT - Unitary Fourier transform
%
% y = ufft(x, ...)
%
% compute the unitary Fourier transform of image. The value of the
% null frequency is equal to mean*sqrt(prod(size(image))). Extra
% argument are passed to fftn function.

  y = fftn(x, varargin{:})/sqrt(prod(size(x)));

end

function x = uifft(x, varargin)
% UIFFT - Unitary inverse Fourier transform
%
% y = uffti(spec, ...)
%
% compute the unitary inverse Fourier transform of spec. The value of
% the null frequency is equal to mean*sqrt(prod(size(spec))). Extra
% argument are passed to ifftn function.

  x = sqrt(prod(size(x)))*ifftn(x, varargin{:});

end

function tf = ir2tf(ir, finalSize)
% IR2TF - Compute the transfert function of the IR
%
% This function make the necessary correct zero-padding, zero
% convention, correct fft2 etc... to compute the transfert function
% of an impulsionnal response.

  sizeir = size(ir);
  irpadded = zeros(finalSize);

  % Zero padding. If you known how to write such things without test
  % the number of dims I'am interessted !
  if (ndims(ir) == 1)
    N = size(ir);
    irpadded(1:N) = ir;
  elseif (ndims(ir) == 2)
    [N M] = size(ir);
    irpadded(1:N, 1:M) = ir;
  elseif (ndims(ir) == 3)
    [N M L] = size(ir);
    irpadded(1:N, 1:M, 1:L) = ir;
  else
    disp('Array of more than 3 dimension are not supported')
    return
  end

  % Zero convention of the fft to avoid the phase problem. Work with
  % odd and even size.
  irpadded = circshift(irpadded, -floor(sizeir/2));
  tf = fftn(irpadded);
end

