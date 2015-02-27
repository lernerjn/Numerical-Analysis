function [ y , tend] = BDF3(h, const2 )
%Author: Jeremy Lerner, Stony Brook University
%Third order implicit BDF solver for the system 
%
%                               [ -k_1   0      0 ]  [y_1]
%      y' = [y_1, y_2 , y_3]' = [  k_1   -k_2   0 ]  [y_2]
%                               [  0      k_2   0 ]  [y_3]
%
%The scheme will stop iterating when either the
% solution does not change within 1e-8 per time step or if the solution
% changes by more than 10 in a single time step (indicating that it
% diverged). Also, this program plots the concentrations from time t=0 to
% time t=tend.
% Input:
%   h: the desired step size
%   const2: k2 in example 9.3, the rate constant for one of the chemical 
%       reactions.
% Output:
%   y: a matrix with all the numerical solutions (at every time step).
%       Where y(1,i) is the concentration of chemical species 1 at time i*h
%   tend: the last time step that the solution was calculated

% rate constant_1 = 1
const1 = 1;
% initial conditions y1(0)=y2(0)=y3(0)=1
y(1:3,1) = 1;

A = [-const1 0 0 ; const1 -const2 0; 0 const2 0];

for i=1:2
    % Since BDF
    % not self starting, we will use a fourth order Runge Kutta method for
    % the first two iterations
    k1 = A*y(1:3,i);
    k2 = A*(y(1:3,i)+(h/2)*k1);
    k3 = A*(y(1:3,i)+(h/2)*k2);
    k4 = A*(y(1:3,i)+h*k3);
    y(1:3,i+1) = y(1:3,i) + (h/6)*(k1+2*k2+2*k3+k4);
end
i = i + 1;
while ( norm(y(1:3,i-1) - y(1:3,i),2) > 1e-8)
    
    
    % BDF, solving the system at each time step (this is an implicit scheme)
     y(1:3,i+1) = ((eye(3,3) -  (6*h/11)*A)) \ ((1/11)*(18*y(1:3,i) - 9*y(1:3,i-1) + 2*y(1:3,i-2))) ;
    
    i = i + 1;
    %check if the solution is diverging
    
    if ( norm(y(1:3,i-1) - y(1:3,i),2) > 10 || i > 1e5)
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
str = sprintf('Plot of the concentrations over time using a third order \n BDF scheme where the rate constant, k2= %g and h = %g', const2, h);
title(str);
str2 = sprintf('BDF3_%g', const2);
legend('species 1','species 2','species 3', 0);
h = figure(1);
saveas(h,str2,'jpg');

end

