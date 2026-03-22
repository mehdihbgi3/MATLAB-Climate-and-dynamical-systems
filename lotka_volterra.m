clear; close all; clc;

t_start = 0;
t_end = 20;
n_steps = 100;
dt = (t_end - t_start) / n_steps;

B0 = 2;
R0 = 1;
y0 = [B0; R0];

time = linspace(t_start, t_end, n_steps + 1);

fprintf(' Lotka-Volterra Predator-Prey Model \n\n');
fprintf('Parameters: t ∈ [%.0f, %.0f], n_steps = %d, dt = %.4f\n', t_start, t_end, n_steps, dt);
fprintf('Initial conditions: B0 = %.1f, R0 = %.1f\n\n', B0, R0);

Y_euler = forward_euler(@lotka_volterra, y0, dt, n_steps);
B_euler = Y_euler(1, :);
R_euler = Y_euler(2, :);

Y_rk4 = rk4(@lotka_volterra, y0, dt, n_steps);
B_rk4 = Y_rk4(1, :);
R_rk4 = Y_rk4(2, :);

figure('Name', 'Time Series Comparison');

subplot(2,2,1);
plot(time, B_euler, 'b-', 'LineWidth', 1.5); hold on;
plot(time, R_euler, 'r-', 'LineWidth', 1.5);
xlabel('Time t');
ylabel('Population');
title(['Forward Euler (dt = ', num2str(dt), ')']);
legend('B(t) - Prey', 'R(t) - Predator', 'Location', 'best');
grid on;

subplot(2,2,2);
plot(time, B_rk4, 'b-', 'LineWidth', 1.5); hold on;
plot(time, R_rk4, 'r-', 'LineWidth', 1.5);
xlabel('Time t');
ylabel('Population');
title(['RK4 (dt = ', num2str(dt), ')']);
legend('B(t) - Prey', 'R(t) - Predator', 'Location', 'best');
grid on;

subplot(2,2,3);
plot(time, B_euler, 'b--', 'LineWidth', 1.5); hold on;
plot(time, B_rk4, 'b-', 'LineWidth', 1.5);
xlabel('Time t');
ylabel('B(t) - Prey');
title('Prey Population Comparison');
legend('Euler', 'RK4', 'Location', 'best');
grid on;

subplot(2,2,4);
plot(time, R_euler, 'r--', 'LineWidth', 1.5); hold on;
plot(time, R_rk4, 'r-', 'LineWidth', 1.5);
xlabel('Time t');
ylabel('R(t) - Predator');
title('Predator Population Comparison');
legend('Euler', 'RK4', 'Location', 'best');
grid on;

sgtitle('Lotka-Volterra: Time Series (n = 100 steps)');

figure('Name', 'Phase Portrait Comparison');

subplot(1,2,1);
plot(B_euler, R_euler, 'b-', 'LineWidth', 1.5); hold on;
plot(B0, R0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(B_euler(end), R_euler(end), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('B (Prey)');
ylabel('R (Predator)');
title(['Forward Euler (dt = ', num2str(dt), ')']);
legend('Trajectory', 'Start', 'End', 'Location', 'best');
grid on;
axis equal;

subplot(1,2,2);
plot(B_rk4, R_rk4, 'b-', 'LineWidth', 1.5); hold on;
plot(B0, R0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(B_rk4(end), R_rk4(end), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('B (Prey)');
ylabel('R (Predator)');
title(['RK4 (dt = ', num2str(dt), ')']);
legend('Trajectory', 'Start', 'End', 'Location', 'best');
grid on;
axis equal;

sgtitle('Phase Portraits (x = B, y = R)');

figure('Name', 'Effect of Time Step Size');

dt_values = [0.2, 0.1, 0.05, 0.01];
colors = {'r', 'b', 'g', 'm'};

subplot(1,2,1);
for i = 1:length(dt_values)
    dt_i = dt_values(i);
    n_i = round((t_end - t_start) / dt_i);
    Y_i = forward_euler(@lotka_volterra, y0, dt_i, n_i);
    plot(Y_i(1,:), Y_i(2,:), colors{i}, 'LineWidth', 1.2); hold on;
end
plot(B0, R0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
xlabel('B (Prey)');
ylabel('R (Predator)');
title('Forward Euler - Different dt');
legend(['dt=', num2str(dt_values(1))], ['dt=', num2str(dt_values(2))], ...
       ['dt=', num2str(dt_values(3))], ['dt=', num2str(dt_values(4))], ...
       'Start', 'Location', 'best');
grid on;

subplot(1,2,2);
for i = 1:length(dt_values)
    dt_i = dt_values(i);
    n_i = round((t_end - t_start) / dt_i);
    Y_i = rk4(@lotka_volterra, y0, dt_i, n_i);
    plot(Y_i(1,:), Y_i(2,:), colors{i}, 'LineWidth', 1.2); hold on;
end
plot(B0, R0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
xlabel('B (Prey)');
ylabel('R (Predator)');
title('RK4 - Different dt');
legend(['dt=', num2str(dt_values(1))], ['dt=', num2str(dt_values(2))], ...
       ['dt=', num2str(dt_values(3))], ['dt=', num2str(dt_values(4))], ...
       'Start', 'Location', 'best');
grid on;

sgtitle('Effect of Time Step Size on Phase Portrait');

figure('Name', 'Effect of Initial Conditions');

IC_values = {[2, 1], [3, 1], [2, 2], [1, 0.5]};
colors_ic = {'b', 'r', 'g', 'm'};

dt_fine = 0.01;
n_fine = round((t_end - t_start) / dt_fine);

for i = 1:length(IC_values)
    y0_i = [IC_values{i}(1); IC_values{i}(2)];
    Y_i = rk4(@lotka_volterra, y0_i, dt_fine, n_fine);
    plot(Y_i(1,:), Y_i(2,:), colors_ic{i}, 'LineWidth', 1.5); hold on;
    plot(y0_i(1), y0_i(2), 'o', 'Color', colors_ic{i}, 'MarkerSize', 8, 'MarkerFaceColor', colors_ic{i});
end
plot(1, 1, 'k*', 'MarkerSize', 15, 'LineWidth', 2);
xlabel('B (Prey)');
ylabel('R (Predator)');
title('Phase Portraits for Different Initial Conditions (RK4)');
legend(['B_0=', num2str(IC_values{1}(1)), ', R_0=', num2str(IC_values{1}(2))], ...
       ['B_0=', num2str(IC_values{2}(1)), ', R_0=', num2str(IC_values{2}(2))], ...
       ['B_0=', num2str(IC_values{3}(1)), ', R_0=', num2str(IC_values{3}(2))], ...
       ['B_0=', num2str(IC_values{4}(1)), ', R_0=', num2str(IC_values{4}(2))], ...
       'Equilibrium (1,1)', 'Location', 'best');
grid on;

figure('Name', 'Conservation Law');

H_euler = B_euler - log(B_euler) + R_euler - log(R_euler);
H_rk4 = B_rk4 - log(B_rk4) + R_rk4 - log(R_rk4);

subplot(1,2,1);
plot(time, H_euler, 'b-', 'LineWidth', 1.5); hold on;
plot(time, H_rk4, 'r-', 'LineWidth', 1.5);
yline(H_euler(1), 'k--', 'LineWidth', 1);
xlabel('Time t');
ylabel('H(B,R)');
title('Conserved Quantity H = B - ln(B) + R - ln(R)');
legend('Euler', 'RK4', 'Initial H', 'Location', 'best');
grid on;

subplot(1,2,2);
plot(time, H_euler - H_euler(1), 'b-', 'LineWidth', 1.5); hold on;
plot(time, H_rk4 - H_rk4(1), 'r-', 'LineWidth', 1.5);
yline(0, 'k--', 'LineWidth', 1);
xlabel('Time t');
ylabel('H(t) - H(0)');
title('Deviation from Initial Value');
legend('Euler', 'RK4', 'Location', 'best');
grid on;

sgtitle('Conservation Law Analysis');

fprintf(' NUMERICAL RESULTS \n\n');
fprintf('Initial H = %.6f\n', H_euler(1));
fprintf('Final H (Euler) = %.6f, Drift = %.6f\n', H_euler(end), H_euler(end) - H_euler(1));
fprintf('Final H (RK4)   = %.6f, Drift = %.6f\n', H_rk4(end), H_rk4(end) - H_rk4(1));



function dydt = lotka_volterra(t, y)
    B = y(1);
    R = y(2);
    
    dydt = zeros(2, 1);
    dydt(1) = B - B * R;
    dydt(2) = -R + B * R;
end

function Y = forward_euler(f, y0, dt, n_steps)
    n_vars = length(y0);
    Y = zeros(n_vars, n_steps + 1);
    Y(:, 1) = y0;
    
    t = 0;
    for n = 1:n_steps
        Y(:, n+1) = Y(:, n) + dt * f(t, Y(:, n));
        t = t + dt;
    end
end

function Y = rk4(f, y0, dt, n_steps)
    n_vars = length(y0);
    Y = zeros(n_vars, n_steps + 1);
    Y(:, 1) = y0;
    
    t = 0;
    for n = 1:n_steps
        k1 = f(t, Y(:, n));
        k2 = f(t + 0.5*dt, Y(:, n) + 0.5*dt*k1);
        k3 = f(t + 0.5*dt, Y(:, n) + 0.5*dt*k2);
        k4 = f(t + dt, Y(:, n) + dt*k3);
        
        Y(:, n+1) = Y(:, n) + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        t = t + dt;
    end
end