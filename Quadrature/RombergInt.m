function [ estimates ] = RombergInt( f,a,b,n )
% Author: Jeremy Lerner, Stony Brook University
% RombergInt estimates the value of the integral of f from a to b using
% Romberg Integration, which uses the Richardson extrapolation based on the
% Trapezoid Rule. Additionally, the program plots the error of the diagonal
% elements of the matrix (those being the best estimates of the function
% value), assuming the true value is pi.
% Input:
%   f: the function to numerically integrate
%   a: the left endpoint
%   b: the right endpoint
%   n: the number of sub intervals to use in the numerical integration
% Output:
%   estimates: N*1 array the best estimates of the integral (on the main
%   diagonal of the matrix). The first element is just trapezoid rule with 
%   one interval, and the ith element uses the Richardson extrapolation on 
%   the previous Richardson extrapolations and so on until it reaches the 
%   trapezoid rule. 

R = zeros(n,n);

%Use the trapezoid rule for the base estimate for every number of
%subintervals
for i = 1:n
    R(i,1) = TrapezoidRule(f,a,b,2^(i-1));    
end

%now use the Richardson extrapolation recursively to refine the estimate
for i =2:n
    for j = 2:i
        R(i,j) = R(i, j-1) + (R(i,j-1) - R(i-1,j-1))/(4^(j-1)-1);
    end
end

estimates = diag(R);
fitRomb = polyfit(log((1:5)'),log( abs(estimates(1:5) - pi)),1);
loglog(1:10, abs(R(1:10)- pi))
xlabel('number of subintervals');
ylabel('log of the error');
str = sprintf('Loglog Plot of the error for Romberg Integration \n Approximate Slope (n=1 to 5) = %g', fitRomb(1));
title(str);
h = figure(1);
%saveas(h,'Romb','jpg')

