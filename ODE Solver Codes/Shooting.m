function [ solns ] = Shooting(  sbelow, sabove)
%Author: Jeremy Lerner, Stony Brook University
%Shooting method to solve the ODE u''=-(1+exp(u)), 0<t<1 . u(0)=0, u(1)=1.
%The program uses bisection to find the initial slope (u'(0) that causes 
%the final value to be 1 (u(1)=1) with successive calls to ode45
%
%Inputs: 
%   sbelow: an initial slope that causes u(1) < 1 , the typical value used
%       in testing was 2
%   sabove: an initial slope that causes u(1) > 1 , the typical value used
%       in testing was 3
%Outputs:
%   solns: a set of column vectors representing the sequential iterations
%   towards the solution. The last column has the final value of u within
%   a tolerance of 1e-8 of 1 (abs(u(1) - 1 )<1e-8)
    
    %The function we want to find a root of is u(end,1) - 1 =0
    fun = @(s)ode45(@sys, [0,1] , [0, s]);
    savg = (sabove+sbelow)/2;
    [~,u] = fun(savg);
    i = 1;
    solns(:,i) = u(:,1);
    
    while ( abs(u(end,1) - 1 )> 1e-3)
        %the value of u(1) is too big, 
        % so throw out the right half of the interval
        if u(end,1) > 1 
            sabove = savg;
            savg =( sabove + sbelow )/2;
            [~,u] = fun(savg);
            i = i + 1;
            solns(:,i) = u(:,1);
        %the value of u(1) is too small,
        % so throw out the left half of the interval
        elseif u(end,1) < 1
            sbelow = savg;
            savg = ( sabove + sbelow )/2;
            [~,u] = fun(savg);
            i = i + 1;
            solns(:,i) = u(:,1);
        %if this is reached, there may be a problem
        else 
            break
        end
    end
    n = size(u,1)-1;
    h = 1/n;
    plot(0:h:1,solns);
    title('Sequential IVP Solutions Until Convergence (to within 1e-3)');
    xlabel('x');
    ylabel('u');

    
end

function [ du ] = sys( t,u )
%The system of first order ODEs that describe the second order ODE:
% u''=-(1+exp(u)). This function basically only exists to be used by ode45

du = zeros(2,1);
du(1)=u(2);
du(2)=-(1+exp(u(1)));

end

