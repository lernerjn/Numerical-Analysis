function [ polyN, polyS, runge ] = CubicSpline_Main( n )
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
%   runge: Runge's function evaluated on the fine mesh in the interval
%   [-1,1]


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
plot((-1:dt:1), polyS, 'g')
hold on
plot((-1:dt:1),runge,'r')
xlabel('t');
ylabel('y');
str = sprintf('Plot of polynomial and cubic spline interpolants for Runge''s function');
title(str);

legend('polynomial', 'cubic spline', 'Runge''s function');

