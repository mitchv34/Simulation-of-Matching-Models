% Define the model
function [dif_eqs, w, profits, vars, params, x_bounds, y_bounds, initial_conditions] = generate_model(assortativity) 
    
    % Define the variables of the system
    syms x y l r % Independent variables
    syms omega_A omega_B sigma_A % Parameters
    syms mu(x) theta(x) % Dependent variables for solution

    % define the part of the production function that depends on x and y
    A = ((omega_A * x^((sigma_A - 1) / sigma_A) + (1 - omega_A) * y^((sigma_A - 1) / sigma_A))^(sigma_A / (sigma_A - 1)));
    % define symbolic expression for th epart that dependes on  l and r
    B = l^omega_B * r^(1-omega_B);

    % Create a symbolic version of the production function
    F = A * B;

    x_bounds = [1.0 100.0]; % Bounds for the support of x Distribution
    y_bounds = [1.0 100.0]; % Bounds for the support of y Distribution

    % Create the distribution objects
    % This is heavier than it shoould be, to allow different distribution
    % of types.
    wuni = makedist('Uniform', 'lower', x_bounds(1), 'upper', x_bounds(2)); 
    funi = makedist('Uniform', 'lower', y_bounds(1), 'upper', y_bounds(2));

    % type_dist = (workers = wuni, firms = funi)
    
    % Create the solver of the model
    H = 100 * (wuni.Upper - wuni.Lower)/(funi.Upper - funi.Lower);

    % Since we are in equilibrium:
    % the variable associated with firms (y) should be substituted with μ
    % l should be θ and r = 1.0 
    f = subs(F, [y , l , r ], [mu, theta, 1.0])	;

    Fx =  subs( diff( F, x ) ,  [y , l , r ], [mu, theta, 1.0]);
    Fxy = subs( diff( F, x, y ) ,  [y , l , r ], [mu, theta, 1.0]);
    Flr = subs( diff( F, l, r ) ,  [y , l , r ], [mu, theta, 1.0]);
    Fxr = subs( diff( F, x, r ) ,  [y , l , r ], [mu, theta, 1.0]);
    Fyl = subs( diff( F, y, l ) ,  [y , l , r ], [mu, theta, 1.0]);

    % Algebraic equations of the model
    w = diff(f, theta(x));
    costs = theta * w; % Costs of the firm
    profits = f - costs; % Profit of the firm
    
    % Creating systems of ordinary differential equations
    if assortativity == "possitive"
        eq1 = diff(theta(x), x) == (H * Fyl - Fxr) / Flr;
        eq2 = diff(mu(x), x) == H / theta(x);
    else
        eq1 = diff(theta(x), x) == -(H * Fyl - Fxr) / Flr;
        eq2 = diff(mu(x), x) == -H / theta(x);
    end
    
    dif_eqs = [eq1, eq2]; % System of differential equations of the model
    
    % Initialize initial conditions
    if  assortativity == "possitive"
        mu_0 = y_bounds(1);
    else 
        mu_0 = y_bounds(2);
    end
    theta_0 = NaN; % Initial condition for θ(x) will be defined later

    initial_conditions = [mu_0, theta_0] ;
    params = [omega_A, omega_B, sigma_A];
    vars = [mu(x), theta(x)];
end