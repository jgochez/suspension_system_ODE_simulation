% Author: Jovanny Gochez
% Date: 5/19/2023
% Description: Linear Damping and Spring Suspension System

clear
clc

% Initialize Constants
const.m = 100; % (kg) mass
const.d = 0.2; % (m) unstretched length of spring
const.k = 4*(10^4); % (N/m) spring coefficient
const.c = 100; % (Ns/m) damping coefficient 
const.g = 9.81; % (m/s^2) gravity

% Initial Conditions
xs_0 = 0.1; % (m) 
dot_xs_0 = 0;% (m/s)
init_cond = [xs_0; dot_xs_0]; % [m, m/s]

% Initialize Simulation Time
t_start = 0; % (s) start
t_end = 20; % (s) end
t_step = 10001; % steps 
t_span = linspace(t_start,t_end,t_step); % timespan 

% Solving our ODE
options = odeset('AbsTol',1e-5,'RelTol',1e-6); % tolerance
[t,y] = ode45(@myODE,t_span,init_cond,options,const);
 
m = const.m;
d = const.d;
k = const.k;
c = const.c;
g = const.g;

% Initialize arrays
Eyw_a = zeros(length(y),1);
fys_a3 = zeros(length(y),1);
fyd_a3 = zeros(length(y),1);

% RECALL: y() comes from @myODE 
for current = 1:length(y)
x_s = y(current,1); %  y = dot_x = [dot_x_1; dot_x_2]
dot_x_s = y(current,2); %  y = dot_x = [dot_x_1; dot_x_2]
Eyw_a(current) = 0.5*m*(dot_x_s^2)+0.5*k*(x_s^2)+(x_s+d)*m*g; 
fys_a3(current) = -k*x_s;
fyd_a3(current) = -c*dot_x_s; 
end

figure % plot solutions

% Total Energy Eyw_a 
subplot(3,2,1)
plot(t,Eyw_a,'Linewidth',3);
hold on
xlabel('Time (s)','fontsize',16,'Interpreter','latex');
ylabel('$E_{yw/a}$ (J)','fontsize',16,'Interpreter','latex');
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',16)
grid on

% Change in Energy delta Eyw_a
subplot(3,2,2)
plot(t,Eyw_a-Eyw_a(1),'Linewidth',3); % Subtract current E with init E
hold on
xlabel('Time (s)','fontsize',16,'Interpreter','latex');
ylabel('$\Delta E_{yw/a}$ (J)','fontsize',16,'Interpreter','latex');
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',16)
grid on

% System force x_s 
subplot(3,2,3)
plot(t,y(:,1),'Linewidth',3); % RECALL: y = dot_x = [dot_x_1; dot_x_2]
hold on
xlabel('Time (s)','fontsize',16,'Interpreter','latex');
ylabel('$x_s$ (m)','fontsize',16,'Interpreter','latex');
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',16)
grid on

% System force dot_x_s
subplot(3,2,4)
plot(t,y(:,2),'Linewidth',3); % RECALL: y = dot_x = [dot_x_1; dot_x_2]
hold on
xlabel('Time (s)','fontsize',16,'Interpreter','latex');
ylabel('$\dot{x}_s$ (m/s)','fontsize',16,'Interpreter','latex');
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',16)
grid on

% Reaction force fyr_b2
subplot(3,2,5)
plot(t,fys_a3,'Linewidth',3);
hold on
xlabel('Time (s)','fontsize',16,'Interpreter','latex');
ylabel('$f^{ys}_{a3}$ (N)','fontsize',16,'Interpreter','latex');
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',16)
grid on

% Reaction force fyr_b3
subplot(3,2,6)
plot(t,fyd_a3,'Linewidth',3);
hold on
xlabel('Time (s)','fontsize',16,'Interpreter','latex');
ylabel('$f^{yd}_{a3}$ (N)','fontsize',16,'Interpreter','latex');
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',16)
grid on


% myODE Function (Must appear at the end of file)
function [dot_x] = myODE(~,x,const)

m = const.m;
c = const.c;
k = const.k;
g = const.g;

% First Order Form
x_1 = x(1); % x_s
x_2 = x(2); % dot_x_s
% x = [x_1; x_2]: for visualization

dot_x_1 = x_2; 
dot_x_2 =((-m*g)/m)-((c*dot_x_1)/m)-((k*x_1)/m) ; % EoM 
dot_x = [dot_x_1; dot_x_2]; 
end
