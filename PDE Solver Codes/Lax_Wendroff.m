function [est, soln,error] = Lax_Wendroff( N, finalT )
% The Lax-Wendroff explicit PDE scheme to solve the 
% partial differential equation v_t + a*u_x = 0 ,
% with the boundary conditions v(x,0) = (sin(pi(x-1)))^2 for 1<=x<=2 and 0
% otherwise and a = 1. Note that the boundary conditions are assumed to be
% periodic for the true and numerical solution. The program uses the 
% explicit Lax-Wendroff scheme for numerically estimating solutions to PDEs
% Input:
%   N: N is the number of intervals in the mesh, meaning that N+1 
%   is the number of points in the mesh.
%   finalT: the time at which to estimate the function v(x,t)
% Output:
%   est: the numerical estimate u(x,finalT)
%   sol: the precalculated analytic solution, this program does not 
%       actually find the true solution, this is hard coded
%   error: a 3*1 array, where error(1) is the 1-norm of the error, error(2)
%   is the 2-norm of the error and error(3) is the infinity norm of the
%   error. Where the error is abs(est - sol).

%finalT = 4;
deltaXSoln = 0.0001;
N1 = 6 / deltaXSoln;


%Setting up the analytic solution

soln = zeros(N1+1, 1);
g = @(x) (sin(pi*(x-finalT-1))).^2;
for k = 1:N1+1
    if ( (mod(k*deltaXSoln - finalT,6) )>= 1 && ( mod(k*deltaXSoln - finalT,6)) <= 2)
        soln(k,1) = g(k*deltaXSoln);
    elseif (  (k*deltaXSoln-finalT) > 2)
        continue;
    end
end

deltaX = 6.0/N;
a = 1;

%The same analytic solution but the same mesh as the numerical solution

solnSmaller = zeros(N+1,1);
for k = 1:N+1
    if ( (k*deltaX - finalT )>= 1 && ( k*deltaX - finalT) <= 2)
        solnSmaller(k,1) = g(k*deltaX);
    elseif (  k*deltaX-finalT > 2)
        continue;
    end
end

%Set deltaT based on the number of mesh points
if (nargin == 3)
    deltaT = deltaTInput;

else
    deltaT = 0.0001;
    if N < 201
        deltaT=0.02;
    elseif N < 401
        deltaT = 0.01;
    elseif N < 801
        deltaT= 0.005;
    elseif N < 1601
        deltaT = 0.0025; 
    end

end


% u(space, time) = u(x,t)
u = zeros(N+1, 1);
u_new = zeros(N+1,1);
f = @(x) (sin(pi*(x-1)))^2;
% set the initial and boundary conditions in the numerical estimate

for k = 1:N+1
    if ( k*deltaX >= 1 && k*deltaX <= 2)
        u(k,1) = f(k*deltaX);
    elseif (  k*deltaX >= 2)
        break;
    end
end

timeSteps = finalT / deltaT;

%Iterate over time
for i = 1:timeSteps 

    %Use periodic boundary conditions to deal with the boundaries
    u_new(1,1) = u(1,1)-(a*deltaT/(2*deltaX))*(u(1+1,1)-u(N,1)) + (((a^2)*(deltaT^2))/(2*deltaX^2))*(u(1+1,1)-2*u(1,1)+u(N,1));
    u_new(N+1,1) = u(N,1)-(a*deltaT/(2*deltaX))*(u(1,1)-u(N-1,1)) + (((a^2)*(deltaT^2))/(2*deltaX^2))*(u(1,1)-2*u(N,1)+u(N+1-1,1));
    
    %Iterate over space, determine each value based on values from the
    %previous time step
    for k = 2:N
        u_new(k,1) = u(k,1)-(a*deltaT/(2*deltaX))*(u(k+1,1)-u(k-1,1)) + (((a^2)*(deltaT^2))/(2*deltaX^2))*(u(k+1,1)-2*u(k,1)+u(k-1,1));
    end
    u = u_new;
end

%output the values for u at the desired finalT
est = u;
%error
error(1) = norm( est - solnSmaller,1);
error(2) = norm( est - solnSmaller,2);
error(3) = norm( est - solnSmaller,inf);

% plot the true solution and the numerical solution
plot((0:deltaX:6), est, '+b')
hold on
plot((0:deltaXSoln:6), soln, 'r')

str = sprintf('Plot of the True and Numerical Solution \n using Lax Wendroff Scheme for N=%f',N);
title(str);
legend('Numerically Calculated Solution', 'True Solution');
xlabel('x');
str2 = sprintf('v(x,%f)',finalT);
ylabel(str2);




end

