
clc 
clear all
close all

% --- Parameters ---
L = 1;            % Length of the rod (spatial domain)
alpha = 1;        % Diffusivity constant
N = 100;          % Number of Fourier terms for analytical solution
t = 0.1;          % Time at which to evaluate solution
n_x = 5;         % Number of spatial points for numerical solution
n_t = 4000;       % Number of time steps for numerical solution
a = 0;            % Start of spatial domain for numerical solution
b = 1;            % End of spatial domain for numerical solution
dx = (b - a) / (n_x + 1); % Spatial step size for numerical solution
dt = t / n_t;     % Time step size for numerical solution (adjusted for final time)
x_analytical = linspace(0, L, N);  % Spatial grid for Fourier series solution
x_numerical = linspace(a, b, n_x + 2)'; % Spatial grid for numerical solution
r = alpha * dt / dx^2;  % Stability parameter for explicit method
k =5;
% --- Initial Condition Function ---
f = @(x) sin(pi*x);  % Example initial condition for both methods

% --- Compute Fourier Coefficients b_n (Analytical Solution) ---
b_n = zeros(1, N);
for n = 1:N
    integrand = @(x) f(x) .* sin(n * pi * x / L);
    b_n(n) = (2 / L) * integral(integrand, 0, L);
end

% --- Compute the Fourier Series Solution (Analytical Solution) ---
u_xt = zeros(size(x_analytical));
for n = 1:N
    u_xt = u_xt + b_n(n) * exp(-((alpha * n * pi / L)^2) * t) .* sin(n * pi * x_analytical / L);
end

% --- Check Stability for Numerical Solution ---
if r > 0.5
    error('The explicit method is unstable, reduce dt or increase n_x.');
end

% --- Numerical Solution (Explicit Finite Difference Method) ---
U_N = zeros(n_x + 2, n_t);  % Solution matrix for numerical solution
U_N(:,1) = f(x_numerical);  % Initial condition for numerical solution
 
U_N(1, :) = 0;              % Boundary condition at x = 0
U_N(end, :) = 0; 
for j = 1:n_t - 1
    for i = 2:n_x+1
        U_N(i,j+1) = r * U_N(i-1,j) + (1 - 2*r) * U_N(i,j) + r * U_N(i+1,j);
    end
end
% Extract numerical solution at t = 0.1
t_step = round(t / dt);  % Find the time step closest to t = 0.1
u_numerical = U_N(:, t_step);

%% --- Interpolation for Error Calculation ---
% Interpolate numerical solution to match analytical grid
u_numerical_interpolated = interp1(x_numerical, u_numerical, x_analytical, 'linear');

% Calculate absolute error
absolute_error = abs(u_numerical_interpolated - u_xt);

% Display the error
fprintf('Absolute Error at t = %.2f:\n', t);
disp(absolute_error);
% --- Plot Both Analytical and Numerical Solutions in One Figure ---
figure;
hold on;

% Plot Analytical Solution (Fourier Series) at t = 0.1
plot(x_analytical, u_xt, 'LineWidth', 2, 'DisplayName', 'Analytical Solution', 'Color', 'b');

% Find the numerical solution closest to t = 0.1
t_step = round(t / dt);  % Find the closest time step to t = 0.1
plot(x_numerical, U_N(:, t_step), 'LineWidth', 2, 'DisplayName', 'Numerical Solution ', 'Color', 'r');

hold off;

% Customize the plot
%title('Analytical vs. Numerical Solution');
xlabel('x');
ylabel('u(x,t)');
legend('show', 'Location', 'best');  % Position the legend at the best location
grid on;
