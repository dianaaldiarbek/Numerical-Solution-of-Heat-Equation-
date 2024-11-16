function [u_xt, U_N, absolute_error] = test_solution(f, L, alpha, N, t, n_x, n_t)
    % Spatial and temporal steps
    dx = L / (n_x + 1);
    dt = t / n_t;
    x_analytical = linspace(0, L, N);
    x_numerical = linspace(0, L, n_x + 2)';
    r = alpha * dt / dx^2;

    % Stability Check
    if r > 0.5
        error('The explicit method is unstable, reduce dt or increase n_x.');
    end

    % --- Compute Fourier Coefficients b_n ---
    b_n = zeros(1, N);
    for n = 1:N
        integrand = @(x) f(x) .* sin(n * pi * x / L);
        b_n(n) = (2 / L) * integral(integrand, 0, L);
    end

    % --- Compute Analytical Solution ---
    u_xt = zeros(size(x_analytical));
    for n = 1:N
        u_xt = u_xt + b_n(n) * exp(-((alpha * n * pi / L)^2) * t) .* sin(n * pi * x_analytical / L);
    end

    % --- Numerical Solution ---
    U_N = zeros(n_x + 2, n_t);
    U_N(:,1) = f(x_numerical);  % Initial condition
    U_N(1, :) = 0;              % Boundary condition at x = 0
    U_N(end, :) = 0;            % Boundary condition at x = L

    for j = 1:n_t - 1
        for i = 2:n_x+1
            U_N(i,j+1) = r * U_N(i-1,j) + (1 - 2*r) * U_N(i,j) + r * U_N(i+1,j);
        end
    end

    % Interpolate and Compute Error
    t_step = round(t / dt);
    u_numerical = U_N(:, t_step);
    u_numerical_interpolated = interp1(x_numerical, u_numerical, x_analytical, 'linear');
    absolute_error = abs(u_numerical_interpolated - u_xt);
end
