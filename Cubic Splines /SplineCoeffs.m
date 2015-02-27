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
