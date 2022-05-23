function acf = computeAutoCorr(y,numlags)
%computeAutoCorr - Compute sample autocorrelation.
%   ACF = computeAutoCorr(Y,NUMLAGS) takes a vector Y of length N, an
%   integer NUMLAGS between 1 and (N-1) and computes sample autocorrelation
%   function of Y. ACF is a vector of length NUMLAGS+1 corresponding to
%   lags 0,1,2,...,NUMLAGS. The first element of ACF is unity (i.e., ACF(1)
%   = 1 at lag 0).

% Reference:
%
%   [1] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
%       Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
%       NJ: Prentice-Hall, 1994.

% Copyright 2016 The MathWorks, Inc.

    % 1. Ensure the sample data is a vector.
    y = y(:);
    N = length(y);

    % 2. numlags must be <= (N-1).
    numlags = min(numlags,N-1);

    % 3. Compute ACF.
    % Convolution, polynomial multiplication, and FIR digital filtering are all
    % the same operation. The FILTER command could be used to compute the ACF
    % (by convolving the de-meaned y with a flipped version of itself), but
    % FFT-based computation is significantly faster for large data sets.
    % The ACF computation is based on [1], pages 30-34, 188:

    nFFT = 2^(nextpow2(length(y))+1);
    
    F = fft(y-mean(y),nFFT);
    F = F.*conj(F);
    
    acf = ifft(F);
    acf = acf(1:(numlags+1)); % Retain non-negative lags
    acf = acf./acf(1);        % Normalize
    acf = real(acf);
end