close all
clear all
clc

% The experimet is to change the value of s_2 to see the effect on the equilibrium
% Define parameters
params = struct("s_1" , 1, ...
                "s_2" , 0, ... % Initialize s_2 to a ridiculously small value (will be updated later)
                "p" , 2/3, ...
                "b" , 0.1,  ...
                "beta", 0.5, ...
                "delta" , 0.2, ...
                "c" , 0.3,  ...
                "r" , 0.05);
                
% Replicating the table in the paper
s_2_vals = 1.20:0.05:1.40;

% Initialize the equilibrium values
s2 = zeros(numel(s_2_vals), 2);
theta = zeros(numel(s_2_vals), 2);
m = zeros(numel(s_2_vals), 2);
u = zeros(numel(s_2_vals), 2);
gamma = zeros(numel(s_2_vals), 2);
phi = zeros(numel(s_2_vals), 2);
w_11 = zeros(numel(s_2_vals), 2);
w_21 = zeros(numel(s_2_vals), 2);
w_22 = zeros(numel(s_2_vals), 2);

for i = 1:numel(s_2_vals)
    params = update_params(params, "s_2", s_2_vals(i));

    % Solve for cs equilibrium
    cs_eq = cs_equilibrium(params);

    % Solve for es equiibrium
    es_eq = es_equilibrium(params);

    % Store the equilibrium values
    s2(i, 1) = params.s_2; % cs equilibrium
    s2(i, 2) = params.s_2; % es equilibrium
    theta(i, 1) = cs_eq.theta; % cs equilibrium
    theta(i, 2) = es_eq.theta; % es equilibrium
    m(i, 1) = cs_eq.m; % cs equilibrium
    m(i, 2) = es_eq.m; % es equilibrium
    u(i, 1) = cs_eq.u; % cs equilibrium
    u(i, 2) = es_eq.u; % es equilibrium
    gamma(i, 1) = cs_eq.gamma; % cs equilibrium
    gamma(i, 2) = es_eq.gamma; % es equilibrium
    phi(i, 1) = cs_eq.phi; % cs equilibrium
    phi(i, 2) = es_eq.phi; % es equilibrium
    w_11(i, 1) = cs_eq.w_11; % cs equilibrium
    w_11(i, 2) = es_eq.w_11; % es equilibrium
    w_21(i, 1) = cs_eq.w_21; % cs equilibrium
    w_21(i, 2) = es_eq.w_21; % es equilibrium
    w_22(i, 1) = cs_eq.w_22; % cs equilibrium
    w_22(i, 2) = es_eq.w_22; % es equilibrium
end

% Create tables
variable_names = {'s₂', 'θ', 'm(θ)', 'u', 'γ', 'ϕ ', ' w₁₁', 'w₂₁', 'w₂₂'};
% Cross-Partial Equilibrium
table_cs = table( s2(:, 1), theta(:, 1), m(:, 1), u(:, 1),...
    gamma(:, 1), phi(:, 1), w_11(:, 1), w_21(:, 1), w_22(:, 1),...
    'VariableNames', variable_names);
% Ex-post Segmetation Equilibrium
table_es = table( s2(:, 2), theta(:, 2), m(:, 2), u(:, 2),...
    gamma(:, 2), phi(:, 2), w_11(:, 2), w_21(:, 2), w_22(:, 2),...
    'VariableNames', variable_names);

% Print the tables
fprintf('Cross-Partial Equilibrium\n');
disp(table_cs);
fprintf('Ex-post Segmetation Equilibrium\n');
disp(table_es);

% Plot the results
% I will evaluate in more points so It will take longer to run
s_2_vals = 1.20:0.01:1.40;
% Initialize the equilibrium values
s2 = zeros(numel(s_2_vals), 2);
theta = zeros(numel(s_2_vals), 2);
m = zeros(numel(s_2_vals), 2);
u = zeros(numel(s_2_vals), 2);
gamma = zeros(numel(s_2_vals), 2);
phi = zeros(numel(s_2_vals), 2);
w_11 = zeros(numel(s_2_vals), 2);
w_21 = zeros(numel(s_2_vals), 2);
w_22 = zeros(numel(s_2_vals), 2);

for i = 1:numel(s_2_vals)
    params = update_params(params, "s_2", s_2_vals(i));

    % Solve for cs equilibrium
    cs_eq = cs_equilibrium(params);

    % Solve for es equiibrium
    es_eq = es_equilibrium(params);

    % Store the equilibrium values
    s2(i, 1) = params.s_2; % cs equilibrium
    s2(i, 2) = params.s_2; % es equilibrium
    theta(i, 1) = cs_eq.theta; % cs equilibrium
    theta(i, 2) = es_eq.theta; % es equilibrium
    m(i, 1) = cs_eq.m; % cs equilibrium
    m(i, 2) = es_eq.m; % es equilibrium
    u(i, 1) = cs_eq.u; % cs equilibrium
    u(i, 2) = es_eq.u; % es equilibrium
    gamma(i, 1) = cs_eq.gamma; % cs equilibrium
    gamma(i, 2) = es_eq.gamma; % es equilibrium
    phi(i, 1) = cs_eq.phi; % cs equilibrium
    phi(i, 2) = es_eq.phi; % es equilibrium
    w_11(i, 1) = cs_eq.w_11; % cs equilibrium
    w_11(i, 2) = es_eq.w_11; % es equilibrium
    w_21(i, 1) = cs_eq.w_21; % cs equilibrium
    w_21(i, 2) = es_eq.w_21; % es equilibrium
    w_22(i, 1) = cs_eq.w_22; % cs equilibrium
    w_22(i, 2) = es_eq.w_22; % es equilibrium
end

% plot the results
% theta
figure;
plot(s2(:, 1), theta(:, 1), '-o', s2(:, 2), theta(:, 2), '-x', Linewidth = 2);
legend("CS Eq", "ES Eq", 'Location','northwest','Orientation','horizontal');
title("Vacancies Realative to Unemployment");

% m
figure;
plot(s2(:, 1), m(:, 1), '-o', s2(:, 2), m(:, 2), '-x', Linewidth = 2);
legend("CS Eq", "ES Eq", 'Location','northwest','Orientation','horizontal');
title("Arrival Rate for Workers");

% u
figure;
plot(s2(:, 1), u(:, 1), '-o', s2(:, 2), u(:, 2), '-x', Linewidth = 2);
legend("CS Eq", "ES Eq", 'Location','northwest','Orientation','horizontal');
title("Unemployment Rate");

% gamma
figure;
hold on
plot(s2(:, 1), gamma(:, 1), '-o', s2(:, 2), gamma(:, 2), '-x', Linewidth = 2);
yline(params.p)% plot p (fraction of low skilled workers) for reference
legend("CS Eq", "ES Eq", "Fraction of Low-Skilled Workers", 'Location','northwest','Orientation','horizontal');
title("Fraction of Unemployed Workers who are Low-Skill");

% phi
figure;
hold on
plot(s2(:, 1), phi(:, 1), '-o', s2(:, 2), phi(:, 2), '-x', Linewidth = 2);
legend("CS Eq", "ES Eq", 'Location','northwest','Orientation','horizontal');
title("Fraction of Low-Skill Vacancies");

% wages
figure;
hold on
plot(s2(:, 1), w_11(:, 1), '-o', s2(:, 2), w_11(:, 2), '-x', Linewidth = 2);
plot(s2(:, 1), w_21(:, 1), '-o', s2(:, 2), w_21(:, 2), '-x', Linewidth = 2);
plot(s2(:, 1), w_22(:, 1), '-o', s2(:, 2), w_22(:, 2), '-x', Linewidth = 2);
legend("w(s_1, s_1) (CS)", "w(s_1, s_1) (ES)", "w(s_2, s_1) (CS)", "w(s_2, s_1) (ES)", "w(s_2, s_2) (CS)", "w(s_2, s_2) (ES)",...
    'Location','northwest','Orientation','horizontal');ß