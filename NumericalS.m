clc 
clear all
close all

%Numerical solution 

% Define parameters
n_x = 10;                  % Number of spatial points
n_t = 400;                   % Number of time steps
a = 0;                     % Start of spatial domain
b = 1;                     % End of spatial domain
t_0 = 0;                   % Initial time
t_n = 0.1;                % Final time
dx = (b - a) / (n_x + 1);  % Spatial step size
dt = (t_n - t_0) / n_t;    % Time step size
x = linspace(a, b, n_x + 2)'; % Spatial grid, including boundary points
alpha = 1;                 % Diffusivity constant
r = alpha * dt / dx^2;     % Stability parameter for explicit method
k =5;
% Check stability for explicit method
if r > 0.5
    error('The scheme is unstable, reduce dt or increase n_x.');
end

U_N = zeros (n_x+2,n_t);
%initial condation
U_N(:,1) =  x.*(b-x);
  % Initial condition for numerical solution
 
U_N(1, :) = 0;              % Boundary condition at x = 0
U_N(end, :) = 0; 
for j = 1: n_t -1
    for i  = 2: n_x +1
        U_N(i,j+1) = r*U_N(i-1,j)+(1-2*r)*U_N(i,j)+r*U_N(i+1,j);
    end
end

% Visualize the solution
figure;
time_steps_to_plot = round(linspace(1, n_t, 10)); % Select 10 time steps to plot
hold on;
for idx = 1:length(time_steps_to_plot)
    plot(x, U_N(:, time_steps_to_plot(idx)), 'DisplayName', sprintf('t = %.3f', (time_steps_to_plot(idx) - 1) * dt));
end
hold off;
%title('Numerical Solution ');
xlabel('x');
ylabel('U(x, t)');
legend('show');
grid on;