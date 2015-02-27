function [ int ] = MonteCarlo( f,a,b,n )
% Author: Jeremy Lerner, Stony Brook University
% MonteCarlo estimates the value of the integral of f from a to b using
% Monte Carlo integration
% Input:
%   f: the function to numerically integrate
%   a: the left endpoint
%   b: the right endpoint
%   n: the number of points to use in the numerical integration
% Output:
%   int: the estimated value of the integral int(f,a,b) (integral of f
%   from a to b)

%the randomly generated points at which to evaluate f
x = (b-a)*rand(n,1) + a;
int = (b-a) * sum( f(x) ) / n;

end

