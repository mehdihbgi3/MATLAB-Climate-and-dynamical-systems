clear; close all; clc;

S0 = 1361;
sigma = 5.67e-8;
epsilon = 0.6;
C = 2.08e8;

Q = S0 / 4;

dt = 1e7;
t_end = 1e9;
n_steps = t_end / dt;
time = (0:n_steps) * dt;
time_years = time / (365.25 * 24 * 3600);

alpha_const = 0.3;

fprintf('Part 1 & 2: Constant Albedo (alpha = %.1f) \n\n', alpha_const);

T0_warm = 320;
T_warm = forward_euler(@(T) rhs_constant_albedo(T, Q, alpha_const, epsilon, sigma, C), ...
                        T0_warm, dt, n_steps);
fprintf('T0 = %d K: Final temperature T_final = %.2f K\n', T0_warm, T_warm(end));

T0_cold = 280;
T_cold = forward_euler(@(T) rhs_constant_albedo(T, Q, alpha_const, epsilon, sigma, C), ...
                        T0_cold, dt, n_steps);
fprintf('T0 = %d K: Final temperature T_final = %.2f K\n', T0_cold, T_cold(end));

T_eq_const = (Q * (1 - alpha_const) / (epsilon * sigma))^0.25;
fprintf('Analytical equilibrium: T_eq = %.2f K\n\n', T_eq_const);

figure('Name', 'Constant Albedo Results');

subplot(2,1,1);
plot(time_years, T_warm, 'r-', 'LineWidth', 1.5); hold on;
plot(time_years, T_cold, 'b-', 'LineWidth', 1.5);
yline(T_eq_const, 'k--', 'LineWidth', 1);
xlabel('Time [years]');
ylabel('Temperature [K]');
title(['Constant Albedo (\alpha = 0.3, \epsilon = ', num2str(epsilon), ')']);
legend(['T_0 = ', num2str(T0_warm), ' K'], ['T_0 = ', num2str(T0_cold), ' K'], ...
       ['T_{eq} = ', num2str(T_eq_const, '%.1f'), ' K'], 'Location', 'best');
grid on;

subplot(2,1,2);
semilogy(time_years, abs(T_warm - T_eq_const), 'r-', 'LineWidth', 1.5); hold on;
semilogy(time_years, abs(T_cold - T_eq_const), 'b-', 'LineWidth', 1.5);
xlabel('Time [years]');
ylabel('|T - T_{eq}| [K]');
title('Convergence to Equilibrium');
legend(['T_0 = ', num2str(T0_warm), ' K'], ['T_0 = ', num2str(T0_cold), ' K'], 'Location', 'best');
grid on;

fprintf(' Part 3: Temperature-Dependent Albedo \n');
fprintf('alpha_p(T) = 0.5 - 0.2 * tanh((T - 265) / 10)\n\n');

T0_ice = 250;
T_ice = forward_euler(@(T) rhs_variable_albedo(T, Q, epsilon, sigma, C), ...
                       T0_ice, dt, n_steps);
fprintf('T0 = %d K: Final temperature T_final = %.2f K\n', T0_ice, T_ice(end));

T0_warm2 = 270;
T_warm2 = forward_euler(@(T) rhs_variable_albedo(T, Q, epsilon, sigma, C), ...
                         T0_warm2, dt, n_steps);
fprintf('T0 = %d K: Final temperature T_final = %.2f K\n', T0_warm2, T_warm2(end));

figure('Name', 'Temperature-Dependent Albedo Results');

subplot(2,1,1);
plot(time_years, T_ice, 'b-', 'LineWidth', 1.5); hold on;
plot(time_years, T_warm2, 'r-', 'LineWidth', 1.5);
xlabel('Time [years]');
ylabel('Temperature [K]');
title('Temperature-Dependent Albedo: \alpha_p(T) = 0.5 - 0.2 tanh((T-265)/10)');
legend(['T_0 = ', num2str(T0_ice), ' K'], ['T_0 = ', num2str(T0_warm2), ' K'], 'Location', 'best');
grid on;

subplot(2,1,2);
alpha_ice = 0.5 - 0.2 * tanh((T_ice - 265) / 10);
alpha_warm2 = 0.5 - 0.2 * tanh((T_warm2 - 265) / 10);
plot(time_years, alpha_ice, 'b-', 'LineWidth', 1.5); hold on;
plot(time_years, alpha_warm2, 'r-', 'LineWidth', 1.5);
xlabel('Time [years]');
ylabel('Albedo \alpha');
title('Evolution of Albedo');
legend(['T_0 = ', num2str(T0_ice), ' K'], ['T_0 = ', num2str(T0_warm2), ' K'], 'Location', 'best');
grid on;

figure('Name', 'Albedo Function and Energy Balance');

subplot(2,2,1);
T_range = 220:0.5:320;
alpha_T = 0.5 - 0.2 * tanh((T_range - 265) / 10);
plot(T_range, alpha_T, 'b-', 'LineWidth', 2);
xlabel('Temperature [K]');
ylabel('Albedo \alpha_p(T)');
title('Temperature-Dependent Albedo');
grid on;
ylim([0.2, 0.8]);

subplot(2,2,2);
incoming = Q * (1 - alpha_T);
outgoing = epsilon * sigma * T_range.^4;
plot(T_range, incoming, 'r-', 'LineWidth', 2); hold on;
plot(T_range, outgoing, 'b-', 'LineWidth', 2);
xlabel('Temperature [K]');
ylabel('Energy flux [W/m^2]');
title('Energy Balance');
legend('Incoming: Q(1-\alpha)', 'Outgoing: \epsilon\sigma T^4', 'Location', 'best');
grid on;

subplot(2,2,3);
net_flux = incoming - outgoing;
plot(T_range, net_flux, 'k-', 'LineWidth', 2); hold on;
yline(0, 'r--', 'LineWidth', 1);
xlabel('Temperature [K]');
ylabel('Net flux [W/m^2]');
title('Net Energy Flux (dT/dt \propto this)');
grid on;

eq_indices = find(diff(sign(net_flux)));
T_equilibria = zeros(size(eq_indices));
for i = 1:length(eq_indices)
    idx = eq_indices(i);
    T_equilibria(i) = T_range(idx) - net_flux(idx) * (T_range(idx+1) - T_range(idx)) / (net_flux(idx+1) - net_flux(idx));
end
fprintf('\nApproximate equilibrium temperatures:\n');
for i = 1:length(T_equilibria)
    if i <= length(eq_indices)
        idx = eq_indices(i);
        if net_flux(idx) > 0 && net_flux(idx+1) < 0
            stability = 'STABLE';
        else
            stability = 'UNSTABLE';
        end
    end
    fprintf('  T_eq(%d) = %.1f K (%s)\n', i, T_equilibria(i), stability);
end

subplot(2,2,4);
dTdt = (1/C) * net_flux;
plot(T_range, dTdt * 1e6, 'b-', 'LineWidth', 2); hold on;
yline(0, 'r--', 'LineWidth', 1);
if ~isempty(T_equilibria)
    plot(T_equilibria, zeros(size(T_equilibria)), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
end
xlabel('Temperature [K]');
ylabel('dT/dt [\times 10^{-6} K/s]');
title('Phase Portrait');
grid on;

figure('Name', 'Multiple Initial Conditions');

T0_values = [240, 250, 260, 265, 270, 280, 290, 300];
colors = jet(length(T0_values));

for i = 1:length(T0_values)
    T_sol = forward_euler(@(T) rhs_variable_albedo(T, Q, epsilon, sigma, C), ...
                          T0_values(i), dt, n_steps);
    plot(time_years, T_sol, 'Color', colors(i,:), 'LineWidth', 1.5); hold on;
end

xlabel('Time [years]');
ylabel('Temperature [K]');
title('Solutions for Various Initial Conditions (Temperature-Dependent Albedo)');
legend_str = arrayfun(@(x) ['T_0 = ', num2str(x), ' K'], T0_values, 'UniformOutput', false);
legend(legend_str, 'Location', 'best');
grid on;



function dTdt = rhs_constant_albedo(T, Q, alpha, epsilon, sigma, C)
    incoming = Q * (1 - alpha);
    outgoing = epsilon * sigma * T^4;
    dTdt = (1/C) * (incoming - outgoing);
end

function dTdt = rhs_variable_albedo(T, Q, epsilon, sigma, C)
    alpha = 0.5 - 0.2 * tanh((T - 265) / 10);
    incoming = Q * (1 - alpha);
    outgoing = epsilon * sigma * T^4;
    dTdt = (1/C) * (incoming - outgoing);
end

function T = forward_euler(rhs, T0, dt, n_steps)
    T = zeros(1, n_steps + 1);
    T(1) = T0;
    for n = 1:n_steps
        T(n+1) = T(n) + dt * rhs(T(n));
    end
end