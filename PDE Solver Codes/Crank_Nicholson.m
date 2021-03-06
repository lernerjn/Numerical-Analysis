function [est, soln, error] = Crank_Nicholson( N, finalT )
%Author: Jeremy Lerner, Stony Brook University
% The Crank-Nicholson implicit PDE scheme to solve
% the partial differential equation v_t + a*u_x = 0 ,
% with the boundary conditions v(x,0) = (sin(pi(x-1)))^2 for 1<=x<=2 and 0
% otherwise and a = 1. Note that the boundary conditions are assumed to be
% periodic for the true and numerical solution. The program uses the 
% implicit Crank-Nicholson scheme for numerically estimating solutions to PDEs
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


%Setting up the analytic solution
deltaXSoln = 0.0001;
N1 = 6 / deltaXSoln;

soln = zeros(N1+1, 1);
g = @(x) (sin(pi*(x-finalT-1))).^2;
for k = 1:N1+1
    if ( (mod(k*deltaXSoln - finalT,6) )>= 1 && ( mod(k*deltaXSoln - finalT,6)) <= 2)
        soln(k,1) = g(k*deltaXSoln);
    elseif (  (k*deltaXSoln-finalT) > 2)
        continue;
    end
end


%The same analytic solution but the same mesh as the numerical solution
deltaX = 6.0/N;
a = 1;

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

if ( timeSteps == 0)
    timeSteps = 1;
end

r = (1/4) * (a * deltaT)/(deltaX);
%Iterate over time
for i = 1:timeSteps

    %The linear system to solve is Ax=Cu=b, where x = u_new
    C = diag(ones(N+1,1),0) + diag(-r*ones(N,1),1) + diag(r*ones(N,1),-1);
    %periodic boundary conditions
    C(1,N+1) = r;
    C(N+1,1) = -r;

    A = diag(ones(N+1,1),0) + diag(r*ones(N,1),1) + diag(-r*ones(N,1),-1);
    %periodic boundary conditions
    A(1,N+1) = -r;
    A(N+1,1) = r;
    b = C * u;
    
    u_new = gmres(A,b);
    u = u_new;
    
end

%output the values for u at the desired finalT
est = u_new;
%error
error(1) = norm( est - solnSmaller,1);
error(2) = norm( est - solnSmaller,2);
error(3) = norm( est - solnSmaller,inf);

% plot the true solution and the numerical solution
plot((0:deltaX:6), est, '+b')
hold on
plot((0:deltaXSoln:6), soln, 'r')

str = sprintf('Plot of the True and Numerical Solution \n using Crank Nicholson Scheme for N=%f',N);
title(str);
legend('Numerically Calculated Solution', 'True Solution');
xlabel('x');
str2 = sprintf('v(x,%f)',finalT);
ylabel(str2);




end

