function p=cdf(pd,x,uflag)
%CDF Cumulative distribution function for piecewise distribution.
%    P=CDF(OBJ,X) returns an array P of values of the cumulative distribution
%    function (CDF) for the piecewise distribution object OBJ, evaluated
%    at the values in the array X.
%
%    P=CDF(OBJ,X,'upper') returns an array P of values of the upper tail
%    probability for the piecewise distribution object OBJ, evaluated
%    at the values in the array X.
%
%    See also PIECEWISEDISTRIBUTION, PIECEWISEDISTRIBUTION/ICDF, PIECEWISEDISTRIBUTION/RANDOM.

%   Copyright 2006-2017 The MathWorks, Inc. 


% Determine the segment that each point occupies
s = segment(pd,x);

% Invoke the appropriate cdf for each segment
p = nan(size(x),'like',x);
for j=1:max(s(:))
    t = (s==j);
    if any(t(:))
        p(t) = pd.distribution(j).cdf(x(t));
    end
end

if nargin>2
    uflag = convertStringsToChars(uflag);
end
if nargin>2 && isequal(uflag,'upper')
    p = 1-p;
end