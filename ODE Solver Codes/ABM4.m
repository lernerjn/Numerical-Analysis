function [ y , tend] = ABM4(h, const2 )
%Author: Jeremy Lerner, Stony Brook University
%Fourth order Adams-Bashforth predictor with an Adams-Moulton corrector
% solver (a predictor-corrector solver) for the system 
%
%                               [ -k_1   0      0 ]  [y_1]
%      y' = [y_1, y_2 , y_3]' = [  k_1   -k_2   0 ]  [y_2]
%                               [  0      k_2   0 ]  [y_3]
%
% The scheme will stop iterating when either the
% solution does not change within 1e-8 per time step or if the solution
% changes by more than 10 in a single time step (indicating that it
% diverged). Also, this program plots the concentrations from time t=0 to
% time t=tend.
% Input:
%   h: the desired step size
%   const2: k2, the rate constant for a chemical reaction.
% Output:
%   y: a matrix with all the numerical solutions (at every time step).
%       Where y(1,i) is the concentration of chemical species 1 at time i*h
%   tend: the last time step that the solution was calculated

% rate constant_1 = 1
const1 = 1;
% initial conditions y1(0)=y2(0)=y3(0)=1
y(1:3,1) = 1;

A = [-const1 0 0 ; const1 -const2 0; 0 const2 0];

for i=1:3
    % Since Adams-Bashforth predictor with an Adams-Moulton corrector is
    % not self starting, we will use a fourth order Runge Kutta method for
    % the first three iterations
    k1 = A*y(1:3,i);
    k2 = A*(y(1:3,i)+(h/2)*k1);
    k3 = A*(y(1:3,i)+(h/2)*k2);
    k4 = A*(y(1:3,i)+h*k3);
    y(1:3,i+1) = y(1:3,i) + (h/6)*(k1+2*k2+2*k3+k4);
end
i = i + 1;
while ( norm(y(1:3,i-1) - y(1:3,i),2) > 1e-8)
    %Adams-Bashforth predictor:
    y(1:3,i+1) = y(1:3,i) + (h/24) * ( 55*A*y(1:3,i) - 59*A*y(1:3,i-1) + 37*A*y(1:3,i-2) - 9*y(1:3,i-3));
    
    %Adams-Moulton corrector
    y(1:3,i+1) = y(1:3,i) + (h/24)*(9*A*y(1:3,i+1) + 19*A*y(1:3,i) -5*A*y(1:3,i-1) + A*y(1:3,i-2));
    
    i = i + 1;
    %check if the solution is diverging
    if ( norm(y(1:3,i-1) - y(1:3,i),2) > 10 || i > 1e7)
        break
    end
    
end


tend = h*(i-1);

plot(0:h:tend, y(1,:), 'b')
hold on
plot(0:h:tend, y(2,:), 'r')
plot(0:h:tend, y(3,:), 'g')

xlabel('time');
ylabel('concentration');
str = sprintf('Plot of the concentrations over time using fourth order \n Adams Bashforth predictor with an Adams-Moulton corrector \n where the rate constant, k2= %g and h = %g', const2, h);
title(str);
str2 = sprintf('ABM4_%g', const2);
legend('species 1','species 2','species 3', 0);
h = figure(1);
%saveas(h,str2,'jpg');

end

