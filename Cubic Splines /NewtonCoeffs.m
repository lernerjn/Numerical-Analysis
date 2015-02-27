function [ coeffs ] = NewtonCoeffs( T, Y )
% NewtonCoeffs calculates the coefficients of a Newton Polynomial
% interpolant for the data T,Y using the recursive divided difference
% formula.
% Input:
%   T: the independent variable values which correspond to the dependent
%       variable values given in Y
%   Y: where Y(i) is the dependent variable values, at T(i)
% Output:
%   coeffs: the coefficients of the Newton Polynomial, where X(1) is
%   the constant term, X(2) is the coefficient of t, etc.
n = size(T,1);
m = size(T,2);

n = max(n,m);

coeffs = zeros(n,1);
for i=1:n
    coeffs(i,1) = (DivDiff(T(1:i), Y(1:i)));
end
