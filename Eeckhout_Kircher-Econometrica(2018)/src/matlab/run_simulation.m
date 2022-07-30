clear all
close all
clc
% Runing the model with three set of parameter

% param_vals = [ω_A, ω_B, σ_A]
param_vals_1 = [0.25, 0.5, 0.9];
param_vals_2 = [0.50, 0.5, 0.9];
param_vals_3 = [0.75, 0.5, 0.9];


[t_1, y_1] = shooting_method(param_vals_1);
[t_2, y_2] = shooting_method(param_vals_2);
[t_3, y_3] = shooting_method(param_vals_3);

% Make plots
fig_mu = figure();
hold on
plot(t_1, y_1(:, 1), LineWidth=1.25)
plot(t_2, y_2(:, 1), LineWidth=1.25)
plot(t_3, y_3(:, 1), LineWidth=1.25)
title("Matching Function (μ(x))")
legend(["ω_A = 0.25", "ω_A = 0.5", "ω_A = 0.75"], 'Location','northwest','NumColumns',1)
hold off

fig_theta = figure();
hold on
plot(t_1, y_1(:, 2), LineWidth=1.25)
plot(t_2, y_2(:, 2), LineWidth=1.25)
plot(t_3, y_3(:, 2), LineWidth=1.25)
title("Firm Size (θ(x))")
ylim([0 1000])
legend(["ω_A = 0.25", "ω_A = 0.5", "ω_A = 0.75"], 'Location','northeast','NumColumns',1)
hold off

fig_w = figure();
hold on
plot(t_1, y_1(:, 3), LineWidth=1.25)
plot(t_2, y_2(:, 3), LineWidth=1.25)
plot(t_3, y_3(:, 3), LineWidth=1.25)
title("Wages (w)")
legend(["ω_A = 0.25", "ω_A = 0.5", "ω_A = 0.75"], 'Location','northwest','NumColumns',1)
hold off

fig_Pi = figure();
hold on
plot(t_1, y_1(:, 4), LineWidth=1.25)
plot(t_2, y_2(:, 4), LineWidth=1.25)
plot(t_3, y_3(:, 4), LineWidth=1.25)
title("Profits (Π)")
legend(["ω_A = 0.25", "ω_A = 0.5", "ω_A = 0.75"], 'Location','northwest','NumColumns',1)
hold off