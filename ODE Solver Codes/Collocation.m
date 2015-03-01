function [ soln ,x ] = Collocation( n )
%Author: Jeremy Lerner, Stony Brook University
%Collocation method (polynomial estimate) to solve the ODE 
%u''=-(1+exp(u)), 0<t<1 . u(0)=0, u(1)=1.
%This program a polynomial of degree n-1, v, to estimate the true solution.
%That is, v must satisfy the boundary conditions and the ODE at the
%interior points
%Inputs:
%   n: number of interior points
%Outputs:
%   soln: the solution, evaluated at n points, 
%       as a column vector (nx1 array), where the first and last points
%       are the boundary conditions
%   x: the coefficients of the monomial functions that compose the
%   estimated solution, an nx1 array.

%dt
h=1/(n);


%We will be solving the equation Ax=b, where x is the vector of
%coefficients of the monomial basis functions up to degree n-1
A = zeros(n,n);

%left BC: when t=0, we just have phi_0=x1 (the constant)
A(1,1)=1;

%right BC: when t=1, 1^j = 1 where j>=0
A(n,:)=1;

%The nodes, the points where the function v's second derivative will be
%matched to the ODE and the RHS will be evaluated
t = linspace(0,1,n);

%Rows 2 to n-1 (inclusive) of A
%Second derivatives of the monomials, which must satisfy the ODE
% phi  = t^(j-1), therefore phi'' = (j-1)*(j-2)*t^(j-3)
for i=2:n-1
    %the second derivative of the first two monomials is zero, so start at
    %phi = t^2, so phi''=2
    for j=3:n
       % phi'' = (j-1)*(j-2)*t^(j-3)
        A(i,j) = (j-1)*(j-2)*t(i)^(j-3);
    end
end

%initial guess, as the book suggests, a straight line connecting the
%boundary conditions, that is u(0)=0 and u(1)=1, so u=t. For use in 
%evaluating the RHS (so the interior points satisfy the ODE)
u = (0:h:1)';

%Monomial basis functions without coefficients, so u = bases*x , where x is
%a column vector of coefficients for the basis functions
bases = zeros(n,n);
%each row is a point t=t_i
for i=1:n
  %each column is the j-1th degree monomial evaluated at t_i, with
  %coefficients 1. Because we will find and multiply by the coefficients 
  %later
  for j=2:n        
    bases(i,j) = t(i)^(j-1);
  end
end

%Vector of coefficients
x = ones(n,1);
xold = zeros(n,1);

%RHS
b = zeros(n,1);

%Iterate until the solution stops significantly changing 
while norm(x-xold,2) > 1e-8
    xold= x;
    %left BC
    b(1) = 0;
    %right BC
    b(n)=1;
    
    %use the estimated/guessed solution for the RHS
    for i = 2:n-1
        b(i) = -(1 + exp(u(i)));
    end
    
    x = A\b;
    %evaluate the solution, because it will be used in the RHS
    u = bases*x;
end

%the function v, the monomial basis functions with the newly found
%coefficients, x. 
soln = bases*x;

n = size(soln,1) - 1;
h = 1/n;

plot(0:h:1,soln);
title('Collocation Method Numerical Solution');
xlabel('x');
ylabel('u');
