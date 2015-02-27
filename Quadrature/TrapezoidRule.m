function [ int ] = TrapezoidRule( f,a,b,n )
% Author: Jeremy Lerner, Stony Brook University
% TrapezoidRule estimates the value of the integral of f from a to b using
% composite Trapezoid rule
% Input:
%   f: the function to numerically integrate
%   a: the left endpoint
%   b: the right endpoint
%   n: the number of sub intervals to use in the numerical integration
% Output:
%   int: the estimated value of the integral int(f,a,b) (integral of f
%   from a to b)

h = (b - a)/n;


x = (a:h:b)';
int = sum ( ( h/2 )* (f(x(1:end-1)) + f(x(2:end))));