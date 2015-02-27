function [ int ] = SimpsonsRule( f,a,b,n )
% Author: Jeremy Lerner, Stony Brook University
% SimpsonsRule estimates the value of the integral of f from a to b using
% composite Simpson's rule (specifically, the 1/3 rule). Composite
% simpson's rule requires an EVEN number of intervals to converge at a
% reasonable rate. However, with an odd number of intervals, it will still
% converge, just very slowly
% Input:
%   f: the function to numerically integrate
%   a: the left endpoint
%   b: the right endpoint
%   n: the number of sub intervals to use in the numerical integration
% Output:
%   int: the estimated value of the integral int(f,a,b) (integral of f
%   from a to b)

h = (b - a)/n;
%check that the number of subintervals is valid
if ( n < 1)
    int = 0;
elseif (n==1)
    int = ((b-a)/6)*(f(a) + 4*f( (a+b)/2) + f(b));
else
    x = (a:h:b)';
    int = (h/3)*(f(a) + f(b) + 2*sum(f( x(3:2:end-2))) + 4*sum(f(x(2:2:end))) );
end


end
