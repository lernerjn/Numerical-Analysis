function [ est , error] = MoLN( n )
%Author: Jeremy Lerner, Stony Brook University
%Using method of lines, this program solvse the heat equation u_t=u_xx on 
%0<=x<=1 and t>0 with initial condition u(0,x)=sin(pi*x) for 0<x<1 and
%Neumann boundary conditions u'(t,0)=u'(t,1)=0. This generates a stiff ODE, 
%so this program uses the matlab function ode15s, which is for stiff ODEs.
%Additionally, the program plots the lines on which the solution was
%found
%inputs:   
%   n: the number of points (including the two boundary points)
%outputs:
%   est: an M*N array, where N=(1/dx)+1 and M is the number of time steps,
%       which is determined by the ODE solver
%   error: the infinity norm of the error (the maximum error) at time t=0.1
%
%the number of points, where the boundary points ARE included, that is
%x_1=x_n = 0
dx = 1/(n-1);

%exact solution u(t,x)=exp(-pi^2(t))sin(pi*x)
exactSolnFun = @(t,x)exp((-pi^2)*(t))*cos(pi*x);
exactSoln = zeros(n,1);
for i = 1:n
    exactSoln(i) = exactSolnFun(0.1,dx*(i-1));
end

%the initial condition
IC = zeros(n,1);
for i=1:n
    IC(i) = cos(pi*(i-1)*dx);
end

%ode15s is used because the ODE is stiff
[t,est] = ode15s(@sys, [0 ,0.1] , IC);

%the value at each x_i at time t=0.1
endTime = est(end,:)';
error = norm(endTime - exactSoln,inf);

x = linspace(0,1,n);
%New mesh for plotting
mesh2 = zeros(n,length(t));
for i=1:n
   for j=1:length(t)
       mesh2(i,j) = x(i);
   end
end

%3D plot
plot3(t, mesh2, est);
grid on;
title('Solution to Heat Equation (with Neumann BCs) Found by Method of Lines');
xlabel('t');
ylabel('x');
zlabel('u(x,t)');

end
function du = sys(t,u)
%set up the matrix to solve the semi-discretized system at each x_i
n = numel(u);
dx = 1/(n-1);
A = diag(ones(n-1,1),1) - diag(2*ones(n,1),0) + diag(ones(n-1,1),-1);
%ghost point for the neumann boundary conditions, where (u_0 - u_2) /2dx =0
A(1,2) = 2;
A(n,n-1)=2;
A = (1/(dx^2))*A;
du = A*u;
end

