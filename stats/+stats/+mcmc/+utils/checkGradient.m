function [isok,err] = checkGradient(fun,x0,stepsize,gradtol)
%checkGradient - Check analytical gradient vs numerical gradient.
%   [ISOK,ERR] = checkGradient(FUN,X0) takes a function handle FUN, a
%   P-by-1 vector X0 and validates that the analytical gradient of FUN at
%   X0 matches the numerical gradient of FUN at X0. Output ISOK is TRUE if
%   analytical gradient is accurate and is FALSE otherwise. Output ERR is a
%   scalar containing the maximum relative error between the elements of
%   the analytical and numerical gradient vectors. Function handle FUN
%   should be callable like this:
%
%       [F,G] = FUN(X0)
%
%   F is the function value at X0 and G is a P-by-1 analytical gradient
%   vector at X0.
%
%   [...] = checkGradient(FUN,X0,STEPSIZE) also accepts a positive scalar
%   STEPSIZE specifying the step size to use in computing the numerical
%   gradient using central differences.
%
%   [...] = checkGradient(FUN,X0,STEPSIZE,GRADTOL) also accepts a positive
%   scalar GRADTOL specifying the relative tolerance for measuring the
%   agreement between the analytical and numerical gradient.

    % 1. Set default values for stepsize and gradtol.
    if nargin < 3
        stepsize = eps(class(x0))^(1/3);
    end

    if nargin < 4
        gradtol = 1e-6;
    end

    % 2. Get the analytical gradient at x0.
    [~,g] = fun(x0);

    % 3. Get the numerical gradient at x0.
    gNum = classreg.learning.fsutils.Solver.getGradient(fun,x0,stepsize);

    % 4. Compute maximum relative error between g and gNum.
    err = max(abs(gNum - g)./max(1.0,abs(g)));

    % 5. Compare err with gradtol.
    if ( err < gradtol )
        isok = true;
    else
        isok = false;
    end
end