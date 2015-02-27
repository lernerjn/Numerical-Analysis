function [ polyN, polyS, runge ] = CubicSpline( n )
%Author: Jeremy Lerner, Stony Brook University
% This program takes in the number of points to use in the
% interpolation for the interval [-1, 1] of Runge's function,
% then solves for the Newton Polynomial and the cubic spline
% polynomials and outputs them on a very fine mesh (with 10,000 points),
% for the purpose of plotting. The program also plots the interpolants and
% Runge's function on the interval [-1 , 1]. Runge's function:
%       f(t) = 1/(1+25t^2)
% Input:
%   n: the number of points to sample Runge's function and use for the
%   polynomial interpolation and the cubic spline interpolation
% Output:
%   polyN: Newton polynomial function evaluation on the fine mesh 
%       in the interval [-1,1]
%   polyS: cubic spline function evaluation on the fine mesh in 
%       the interval [-1,1]
%   runge: Runge's function evaluated on the fine mesh in the interval [-1,1]


f = @(x) 1./(1+25*x.^2);

T = linspace(-1,1,n)';
Y = zeros(n,1);
Y(1:size(T,1)) = f(T(1:size(T,1)));

X = NewtonCoeffs( T, Y );

N = 10000;
dt = 2/N;
polyN = zeros(N,1);
runge = zeros(N,1);
j=1;
for i =-1:dt:1
    polyN(j) = NewtonPolynomial( T,Y, i, X );
    runge(j) = f(i);
    j=j+1;
end

x = SplineCoeffs( T,Y );
j=1;
k=1;
polyS = zeros(N,1);
index = 1;
for i = -1:dt:1
    
    %check which spline to use in the evaluation
    if ( i >= T(k+1))
        k = k+1;
        index = index+4;
        if index > numel(x)
            break
        end
    end
    %g is the spline for the given interval
    g = @(t)x(index)+x(index+1)*t + x(index+2)*t.^2 + x(index+3)*t.^3;
    polyS(j) = g(i);
    j = j+1;

end

polyS(N+1,1) = g(1);
plot((-1:dt:1),polyN)
hold on
plot((-1:dt:1), polyS, 'g')

plot((-1:dt:1),runge,'r')
xlabel('t');
ylabel('y');
str = sprintf('Plot of polynomial and cubic spline interpolants for Runge''s function');
title(str);

legend('polynomial', 'cubic spline', 'Runge''s function');

end

function [ x] = SplineCoeffs( T,Y )
% SplineCoeffs calculates the coefficients for a cubic spline interpolation
% for the independent variable T and dependent variable Y. 
% Input:
%   T: the independent variable values which correspond to the dependent
%       variable values given in Y
%   Y: where Y(i) is the dependent variable values, at T(i)
% Output:
%   x: the coefficients of a cubic spline that interpolates the given data.
%   Where p = x(i) + t*x(i+1) + t^2*x(i+2) + t^3*x(i+3), for i = 1,5,9,...
%   These are the coefficients for the spline in the given intervals, so
%   i=1 corresponds to the interval T(1) to T(2), and i=2 corresponds to
%   the interval T(2) to T(3), etc.

n = max(size(T,1), size(T,2));

%the constant terms in the function, for left and right
A = zeros(4*(n-1));

%function value on the left
row = 1;
index = 1;
for i =1:4:4*(n-1)
    A(row,i) = 1;
    A(row,i+1) = T(index);
    A(row,i+2) = T(index)^2;
    A(row,i+3) = T(index)^3;
    index = index + 1;
    row = row + 2;

end

%function value on the right
index = 1;
row = 2;
for i =1:4:4*(n-1)
    A(row,i)=1;
    A(row,i+1) = T(index+1);
    A(row,i+2) = T(index+1)^2;
    A(row,i+3) = T(index+1)^3;
    index = index + 1;
    row = row + 2;

end

%first derivatives
index = 2;
row = 2*(n-1)+1;
for i=1:4:4*(n-1)-4
    A(row,i+1) = 1;
    A(row,i+2) = 2*T(index);
    A(row,i+3) = 3*T(index)^2;
    row = row + 1;
    index = index + 1;

end

%first derivatives
index = 2;
row = 2*(n-1)+1;
for i=5:4:4*(n-1)
    A(row,i+1) = -1;
    A(row,i+2) = -2*T(index);
    A(row, i+3) = -3*T(index)^2;
    row = row + 1;
    index = index + 1;
end

%second derivatives
index = 2;
row = 3*(n-1);
for i=1:4:4*(n-1)
    A(row,i+2) = 2;
    A(row,i+3) = 6*T(index);
    row = row + 1;
    index = index + 1;
end

%second derivatives
index = 2;
row = 3*(n-1);
for i=5:4:4*(n-1)
    A(row,i+2) = -2;
    A(row,i+3) = -6*T(index);
    row = row + 1;
    index = index + 1;
end

A(end - 1, : ) = 0;
A(end, :) =0;
%second derivatives zero at the endpoints (natural spline)
A(end - 1, 3) = 2;
A(end-1,4) = 6*T(1);

A(end, end-1) = 2;
A(end,end) = 6*T(end);

%The right hand side. Where Ax=b
b = zeros(4*(n-1),1);
loopVar = 1;
b(1) = Y(1);
for i =1:2:4*(n-1)/2-2
    b(i+1)=Y(loopVar+1);
    b(i+2)=Y(loopVar+1);
    loopVar = loopVar + 1;
end
b(4*(n-1)/2) = Y(end);

x = A\b;
end


function [ out ] = NewtonPolynomial( T,Y, t, X )
%Author: Jeremy Lerner, Stony Brook University
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
end

function [ coeffs ] = NewtonCoeffs( T, Y )
%Author: Jeremy Lerner, Stony Brook University
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
end

function [ x ] = DivDiff( T, Y )
%Author: Jeremy Lerner, Stony Brook University
% Divided differences (using the recursive formula),
% used to find the Newton Polynomial
if ( size(T) == 1)
    x = Y;
else
    x = (DivDiff(T(2:end), Y(2:end)) - DivDiff(T(1:end-1),Y(1:end-1)))/(T(end) - T(1));
end
end


