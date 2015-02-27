function [ out ] = NewtonPolynomial( T,Y, t, X )
% NewtonPolynomial calculates the function value of a Newton Polynomial
% interpolant for the data T,Y at the location t. An optional input is a
% vector of coefficients for a Newton Polynomial.
% for the independent variable T and dependent variable Y. 
% Input:
%   T: the independent variable values which correspond to the dependent
%       variable values given in Y
%   Y: where Y(i) is the dependent variable values, at T(i)
%   t: the point to find the value of the Newton Polynomial
%   X: OPTIONAL: the coefficients of the Newton Polynomial, where X(1) is
%   the constant term, X(2) is the coefficient of t, etc.
% Output:
%   out: the function evaluation at the point t, that is p(t)

if (nargin == 3)
    X = NewtonCoeffs(T,Y);
end

n = size(T,1);
m = size(T,2);

n = max(n,m);


out = X(1);
for i=2:n
    mid = 1;
    for j = 1:i-1
        mid = mid*(t-T(j));
    end
    out = out + X(i)*mid;
end

