%% 
close all
clear all
clc


% Define parameters
L = 1;            % Length of the rod (spatial domain)
alpha = 1;        % Diffusivity constant
N = 100;           % Number of Fourier terms
t = 0.1;          % Time at which to evaluate solution
x = linspace(0, L, N);  % Spatial grid
k =5;
% Define initial condition as a function handle
f = @(x) x.*(L-x);  % Example initial condition

% Compute Fourier coefficients b_n
b_n = zeros(1, N);
for n = 1:N
    integrand = @(x) f(x) .* sin(n * pi * x / L);
    b_n(n) = (2 / L) * integral(integrand, 0, L);
end

% Evaluate the Fourier series solution u(x, t)
u_xt = zeros(size(x));
for n = 1:N
    u_xt = u_xt + b_n(n) * exp(-(alpha* n * pi / L)^2 * t) .* sin(n * pi * x / L);
end

% Plot the solution
plot(x, u_xt, 'LineWidth', 2);
xlabel('x');
ylabel('u(x, t)');
title(['Solution of the heat equation at t = ' num2str(t)]);
grid on;

disp(b_n)