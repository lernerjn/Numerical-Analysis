function [ soln ] = FiniteDiff( n )
%Author: Jeremy Lerner, Stony Brook University
%Finite Difference method to solve the ODE 
%u''=-(1+exp(u)), 0<t<1 . u(0)=0, u(1)=1.
%This program uses a mesh with n nodes to discretize and solve the above
%ODE. The Matlab nonlinear equation solver, fsolve, is used for the system
%of equations that is generated
%Inputs:
%   n: number of interior points
%Outputs:
%   soln: the solution, evaluated at n+2 points, 
%       as a column vector (n+2x1 array), where the first and last points
%       are the boundary conditions

%dt
h=1/(n+1);

%initial guess, as the book suggests, a straight line connecting the
%boundary conditions, that is u(0)=0 and u(1)=1, so u=t
x0 = (h:h:(1-h))';

%The set of nonlinear equations to solve, with a specific n and h, as
%defined above. The function is defined in this same file, below
fun = @(u)FinDiffEqns(u,n,h);
%Call the built in nonlinear equation solver
[fval] = fsolve(fun,x0);
%since n represents the number of interior points, and we solve for only
%the interior points in the above line, add the boundary values to the
%solution vector
soln = [0; fval ; 1];

plot(0:h:1,soln);
title('Finite Difference Numerical Solution');
xlabel('x');
ylabel('u');

end


function [ F ] = FinDiffEqns( u,n,h)
    %The system of nonlinear equations for the interior points of the ODE above
    F = zeros(n,1);

    %Just one interior point, so we need to use the boundary conditions for
    %both sides
    if (n==1)
        F(1)= (1 - 2*u(1) + 0)/(h^2) + 1 + exp(u(1));
    else
        %left boundary condition is involved in the formula, u_0=u(0)=0
        F(1)= (u(2) - 2*u(1) + 0)/(h^2) + 1 + exp(u(1));
        for i=2:n-1
            F(i) = (u(i+1) - 2*u(i) + u(i-1))/(h^2) + 1 + exp(u(i));
        end
        %right boundary condition is involved in the formula, u_(n+1)=u(1)=1
        F(n)=(1 - 2*u(n) + u(n-1))/(h^2) + 1 + exp(u(n));
    end

end
