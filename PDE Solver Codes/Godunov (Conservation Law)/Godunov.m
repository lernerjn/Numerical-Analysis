function [x, u] = Godunov( dtInv )
%Author: Jeremy Lerner, Stony Brook University
% Godunov's method for numerical conservation laws to plot the solution to 
% Burger's equation where u_t + u*u_x = 0 , t>0 at t=2, with initial condition
%           { -1   , x < 1
% u_(x,0)=  {
%           { 1    , x > 1
%
% We will use the numerical conservation method: 
%     U(j,n+1) = U(j,n) - (dt/dx)( F(U(j,n), U(j+1,n) - F(U(j,n), U(j+1,n))
%
%           { f(v)   , (f(v)-f(w))/(v-w) >= 0
% F(v,w) =  {
%           { f(w)   , (f(v)-f(w))/(v-w) < 0
%
%   dt/dx = 1/2 -> dt = dx/2
%
% inputs:
%       dtInv: the inverse of the time step, 1/dt
% outputs:
%       x: the mesh from 0 to 2 
%       u: the solution evaluated at all the points in x


% Note: this method converges for dtInv = 2n but does not converge for
% dtInv = 2n+1 (when dtInv = 2n+1, an odd number, the numerical solution
% gives a standing wave), which is a property of this conservation law method.


dt = 1/dtInv;
%left end point
a = 0;
%right end point
b = 2;

%beginning time, t=0
%end time
t_end = 2;

%dt/dx = 1/2
dx = 2*dt;
%dx = 1/(n-1), where n is the total number of points (inlcuding boundary
%points)

x = [a:dx:b]';
n = numel(x);


%cell averaged initial condition
u = ones(n,1);
for i=1:n
    if ( x(i) < 1 )
        %integral of the constant function u=1 from x(i) to x(i+1) divided by dx
        u(i,1) = -1 ;
    elseif (x(i) > 1)
        u(i,1) = 1 ;
    end
end

%if there is an even number of points, the cell boundary for the left and
%right end at the discontinuity and will be entirely determined by the left
%and right side respectively. However, if there is an odd number of points,
%then the middle point will be equally influenced by the left and right
if (mod(n,2)== 1)
    u(ceil(n/2)) = 0;
end
unew = zeros(n,1);

%go to the end time of t_end
for i=0:dt:t_end
    
    for j = 2: n-1
        unew(j , 1) = u(j,1) - (1/2)*(Flux(u(j,1),u(j+1)) - Flux(u(j-1,1),u(j,1)));
    end
    u = unew;
    u(1,1) = -1;
    u(n,1) = 1;
    
end

plot(x,u);
title(sprintf('Godunov''s method (a numerical conservation law) used to solve Burger''s equation,\nSolution shown at t=2'));
xlabel('x');
ylabel('u');

end



%numerical flux, using the Godunov method:
function F = Flux(v, w)
    f = @(u) (1/2)*u.^2;
    if v == w
        F = f(v);
    elseif ((f(v) - f(w))/(v-w)) >=0
        F = f(v);
    else
        F = f(w);
    end
    
    
end    

