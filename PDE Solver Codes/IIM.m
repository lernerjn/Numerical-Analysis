function [ u, error , soln] = IIM( N )
%Author: Jeremy Lerner, Stony Brook University
%Immersed Interface Method for uxx = f(x), u(-1)=-1/e 
% Randall Leveque's Immersed Interface Method to solve the following
% elliptic problem:
% u_{xx}=f(x)
%  , u(-1)=\frac{-1}{e}
%  , u(1)=0
% f(x)=-6x*e^{-x^{2}}+4x^{3}*e^{-x^{2}}
%   , x\in[-1,0)
% f(x)=-6x
%   , x\in(0,1]
%  
% u\left(0^{+}\right)-u\left(0^{-}\right)=1
%  , u_{x}\left(0^{+}\right)-u_{x}\left(0^{-}\right)=-1
% 
% The exact solution is 
% u(x)=x*e^{-x^{2}}
%   , x\in[-1,0)
%  
% u(x)=1-x^{3}
%   , x\in(0,1]
% Inputs: 
%   N: Number of points in the grid, N must be even to assure that x=0 is
%   in the middle of a computational cell, as desired. If N is odd, N+1
%   will be used as the number of grid points
% Outputs:
%   u: numerical solution (using N evenly spaced points -1 to 1)
%   error: the error vector (the error at each coordinate), note this is
%      the error, not the absolute value of the error, that is numerical
%       solution minus true solution
%   soln: the true solution (using N evenly spaced points -1 to 1)

%Force N to be even
if (mod(N,2) ~= 0 )
    N = N+1;
end

%Set up the grid for the x coordinates, since there are N points, there are
%N-1 subintervals, meaning h=(b-a)/(N-1)
%a = x(left) = -1 , b = x(right) = 1
h = (1 - (-1))/(N-1);
grid = linspace(-1,1,N);



%Setting up Au=b
A = diag(1*ones(1,N-1),-1) + diag(-2*ones(1,N),0) + diag(1*ones(1,N-1),1);
%Left boundary condition: u1 = (1/h^2) * (u0 - 2u1 +u2), u0 is the first
%point and u0 = -1/e
A(1,:) = 0;
A(1,1) = 1;
%Right boundary condition, u_N = 0
A(N,:) = 0;
A(N,N) = 1;
A = A*(1/h^2);

%right hand side
b = zeros(N,1);
%left side boundary condition
%b(1) =  -6*grid(1)*exp(-(grid(1).^2)) + 4*(grid(1).^3)*exp(-(grid(1).^2)) - (-exp(-1))/h^2;
b(1) = (-exp(-1))/h^2 ;
for i = 2:N/2
    b(i) = -6*grid(i)*exp(-(grid(i).^2)) + 4*(grid(i).^3)*exp(-(grid(i).^2));
end

for i = ((N/2)+1):(N-1)
    b(i) = -6*grid(i);
end
%right side boundary condition
b(N) = 0;
%Adding C_j and C_j+1 to the RHS


%                       C_j
%The folllwing line is true, but the zeros in the last term are
%unnecessary for the math, but helpful to explain the thought process
%b(N/2) = b(N/2) + (1/h^2)*1*(1+grid(N/2+1)*(-1) - ((grid(N/2+1)^2)/2 )* (0 + 0));
%However, this line has no unnecessary function calls or calculations
b(N/2) = b(N/2) + (1/h^2)*1*(1+grid(N/2+1)*(-1));


%                       C_j+1
%The folllwing line is true, but the zeros in the last term are
%unnecessary for the math, but helpful to explain the thought process
%b(N/2 + 1) = b(N/2+1) + (1/h^2)*1*(-1 - grid(N/2)*(-1) - (((grid(N/2))^2)/2)*(0 - 0 + 0));
%However, this line has no unnecessary function calls or calculations
b(N/2 + 1) = b(N/2 + 1) + (1/h^2)*1*(-1 - grid(N/2)*(-1));

%Now, around x = 0:
%system for coeffs gammas is: A= [1 1 1; x(j-1) x(j) x(j+1) ;x(j-1)^2
%x(j)^2 x(j+1)^2 ] * [gamma1 gamma2 gamma3]' = [0 0 2]. Thus 
%gamma1 = 1/h^2; gamma2 = -2/h^2 and gamma3 = 1/h^2, therefore the system
%around 0 is the same as the system on the whole
u =A\b;

%Plot the numerical solution
plot(linspace(-1,1,N),u)
hold on

x = linspace(-1,1,N);
soln = zeros(N,1);
%Piece together the true solution
for i = 1:N
    
    if ( x(i) < 0 )
        %To the left of zero
        soln(i) = x(i) * exp(-(x(i)^2));
    
    else
        %To the right of zero
        soln(i) = 1 - x(i)^3;
    end
    
end


error = u - soln;
plot( x, soln, 'or');

plot(x, error, 'og')

xlabel('x');
ylabel('u(x)');
str = sprintf('Plot of the numerical solution found \n using Leveque''s Immersed Interface Method with N=%g' ,N);
title(str);
h = figure(1);
str2 = sprintf('IIM%g', N);
legend('Numerical solution' , 'Analytic solution' , 'error', 0);
% saveas(h,str2,'jpg');

end

