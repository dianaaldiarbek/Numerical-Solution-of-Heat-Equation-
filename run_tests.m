% --- Parameters ---
L = 1; alpha = 1; N = 100; t = 0.1; n_x = 15; n_t = 4000;

% Define initial conditions
f_list = {
    @(x) sin(pi*x),       % Test 1: Sine function
    @(x) x .* (1-x),      % Test 2: Parabolic function
    @(x) exp(-5*(x-0.5).^2), % Test 3: Gaussian function
};

test_names = {'Sine Function', 'Parabolic Function', 'Gaussian Function'};

% Run Tests
for k = 1:length(f_list)
    fprintf('------------------------------------\n');
    fprintf('Test %d is starting: %s\n', k, test_names{k});
    tic;  % Start timing
    
    % Extract test function
    f = f_list{k};
    [u_xt, U_N, absolute_error] = test_solution(f, L, alpha, N, t, n_x, n_t);
    elapsed_time = toc;  % Stop timing

    % Display results
    max_error = max(absolute_error);
    fprintf('CPU time = %.6f sec\n', elapsed_time);
    fprintf('Maximum Absolute Error = %.6e\n', max_error);

    % Check if test passed
    if max_error < 1e-3  % Define a tolerance for passing
        fprintf('Test %d passed successfully !!!\n', k);
    else
        fprintf('Test %d FAILED !!!\n', k);
    end
end
fprintf('------------------------------------\n');
