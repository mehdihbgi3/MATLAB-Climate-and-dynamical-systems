clear; close all; clc;

sigma = 10;
r = 28;
b = 8/3;
params = [sigma, r, b];

dt = 0.01;
n_steps = 16000;
n_loc = 10;
n_assim = n_steps / n_loc;

x0_true = [8; 0; 30];
per = [5; 5; 5];
x0_forecast = x0_true + per;

obs_var = 2.0;
R = obs_var * eye(3);

forecast_var = 8.0;
Pf = forecast_var * eye(3);

H = eye(3);

time = (0:n_steps) * dt;
time_assim = (0:n_loc:n_steps) * dt;

X_truth = zeros(3, n_steps + 1);
X_truth(:, 1) = x0_true;

X_forecast = zeros(3, n_steps + 1);
X_forecast(:, 1) = x0_forecast;

X_analysis = zeros(3, n_steps + 1);
X_analysis(:, 1) = x0_forecast;

Y_obs = zeros(3, n_assim + 1);

for k = 1:n_steps
    X_truth(:, k+1) = rk4(@lorenz, time(k), X_truth(:, k), dt, params);
end

figure('Name', 'Lorenz Attractor');

subplot(1, 2, 1);
plot3(X_truth(1,:), X_truth(2,:), X_truth(3,:), 'b-', 'LineWidth', 0.5);
xlabel('x'); ylabel('y'); zlabel('z');
title('Lorenz Attractor - 3D View');
grid on;
view(45, 30);

subplot(1, 2, 2);
plot(X_truth(1,:), X_truth(3,:), 'b-', 'LineWidth', 0.5);
xlabel('x'); ylabel('z');
title('Lorenz Attractor - xz Projection');
grid on;

x_current = x0_forecast;

assim_idx = 1;
for k = 1:n_steps
    current_time = time(k);
    
    if mod(k-1, n_loc) == 0
        Y_obs(:, assim_idx) = generate_observations(X_truth(:, k), R);
        
        x_analysis = OI(x_current, Y_obs(:, assim_idx), H, Pf, R);
        
        X_analysis(:, k) = x_analysis;
        X_forecast(:, k) = x_current;
        
        x_current = x_analysis;
        
        assim_idx = assim_idx + 1;
    else
        X_analysis(:, k) = x_current;
        X_forecast(:, k) = x_current;
    end
    
    x_current = rk4(@lorenz, current_time, x_current, dt, params);
end

Y_obs(:, end) = generate_observations(X_truth(:, end), R);
x_analysis = OI(x_current, Y_obs(:, end), H, Pf, R);
X_analysis(:, end) = x_analysis;
X_forecast(:, end) = x_current;

t_start = 120;
t_end = 160;
idx_range = find(time >= t_start & time <= t_end);

obs_time_indices = 1:n_loc:n_steps+1;
obs_in_window = (time_assim >= t_start) & (time_assim <= t_end);
obs_times_plot = time_assim(obs_in_window);
obs_indices_plot = obs_time_indices(obs_in_window);

figure('Name', 'Data Assimilation Results - z Component');

subplot(3, 1, 1);
plot(time(idx_range), X_truth(3, idx_range), 'r-', 'LineWidth', 1.5); hold on;
plot(time(idx_range), X_analysis(3, idx_range), 'b-', 'LineWidth', 1);
plot(obs_times_plot, Y_obs(3, obs_in_window), 'k*', 'MarkerSize', 5);
legend('Truth', 'Analysis', 'Observations', 'Location', 'best');
ylabel('z');
title('Truth vs. Analysis vs. Observations');
xlim([t_start, t_end]);
grid on;

subplot(3, 1, 2);
plot(time(idx_range), X_forecast(3, idx_range) - X_truth(3, idx_range), 'b-', 'LineWidth', 1);
ylabel('Error');
title('Forecast - Truth');
xlim([t_start, t_end]);
grid on;

subplot(3, 1, 3);
plot(time(idx_range), X_analysis(3, idx_range) - X_truth(3, idx_range), 'b-', 'LineWidth', 1);
xlabel('Time');
ylabel('Error');
title('Analysis - Truth');
xlim([t_start, t_end]);
grid on;

sgtitle('Figure 1: Error for the z-component of the Lorenz model');

figure('Name', 'Data Assimilation Results - All Components');
comp_labels = {'x', 'y', 'z'};

for comp = 1:3
    subplot(3, 3, (comp-1)*3 + 1);
    plot(time(idx_range), X_truth(comp, idx_range), 'r-', 'LineWidth', 1.5); hold on;
    plot(time(idx_range), X_analysis(comp, idx_range), 'b-', 'LineWidth', 1);
    plot(obs_times_plot, Y_obs(comp, obs_in_window), 'k*', 'MarkerSize', 4);
    ylabel(comp_labels{comp});
    xlim([t_start, t_end]);
    if comp == 1
        title('Truth vs Analysis');
        legend('Truth', 'Analysis', 'Obs', 'Location', 'best');
    end
    grid on;
    
    subplot(3, 3, (comp-1)*3 + 2);
    plot(time(idx_range), X_forecast(comp, idx_range) - X_truth(comp, idx_range), 'b-', 'LineWidth', 1);
    xlim([t_start, t_end]);
    if comp == 1, title('Forecast - Truth'); end
    grid on;
    
    subplot(3, 3, (comp-1)*3 + 3);
    plot(time(idx_range), X_analysis(comp, idx_range) - X_truth(comp, idx_range), 'b-', 'LineWidth', 1);
    xlim([t_start, t_end]);
    if comp == 1, title('Analysis - Truth'); end
    if comp == 3, xlabel('Time'); end
    grid on;
end

rms_forecast = sqrt(mean((X_forecast - X_truth).^2, 2));
rms_analysis = sqrt(mean((X_analysis - X_truth).^2, 2));

fprintf('RMS Errors:\n');
fprintf('Component    Forecast    Analysis    Improvement\n');
fprintf('-------------------------------------------------\n');
for i = 1:3
    improvement = (1 - rms_analysis(i)/rms_forecast(i)) * 100;
    fprintf('    %s        %7.4f     %7.4f       %5.1f%%\n', ...
        comp_labels{i}, rms_forecast(i), rms_analysis(i), improvement);
end

fprintf('\nMean Absolute Errors:\n');
mae_forecast = mean(abs(X_forecast - X_truth), 2);
mae_analysis = mean(abs(X_analysis - X_truth), 2);
fprintf('Component    Forecast    Analysis\n');
fprintf('----------------------------------\n');
for i = 1:3
    fprintf('    %s        %7.4f     %7.4f\n', ...
        comp_labels{i}, mae_forecast(i), mae_analysis(i));
end

fprintf('\nCorrelation with Truth:\n');
fprintf('Component    Forecast    Analysis\n');
fprintf('----------------------------------\n');
for i = 1:3
    corr_f = corrcoef(X_forecast(i,:), X_truth(i,:));
    corr_a = corrcoef(X_analysis(i,:), X_truth(i,:));
    fprintf('    %s        %7.4f     %7.4f\n', ...
        comp_labels{i}, corr_f(1,2), corr_a(1,2));
end

function dxdt = lorenz(params, state)
    sigma = params(1);
    r     = params(2);
    b     = params(3);
    
    x = state(1);
    y = state(2);
    z = state(3);
    
    dxdt = zeros(3, 1);
    dxdt(1) = sigma * (y - x);
    dxdt(2) = r * x - y - x * z;
    dxdt(3) = x * y - b * z;
end

function x_next = rk4(f, t_n, x_n, dt, params)
    k1 = f(params, x_n);
    k2 = f(params, x_n + 0.5 * dt * k1);
    k3 = f(params, x_n + 0.5 * dt * k2);
    k4 = f(params, x_n + dt * k3);
    
    x_next = x_n + dt * (k1 + 2*k2 + 2*k3 + k4) / 6;
end

function x_a = OI(x_f, y, H, Pf, R)
    innovation = y - H * x_f;
    S = H * Pf * H' + R;
    K = (Pf * H') / S;
    x_a = x_f + K * innovation;
end

function y = generate_observations(x_true, R)
    n = length(x_true);
    L = chol(R, 'lower');
    noise = L * randn(n, 1);
    y = x_true + noise;
end