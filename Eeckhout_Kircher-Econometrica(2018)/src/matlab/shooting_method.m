% This function iteratively update the initial guess until convergence
% This is a shooting method to solve the IVP

function [t,y] = shooting_method(param_vals)

    % Set tolerance and maximum number of iterations
    tol = 1e-5;
    max_iter = 100;

    % Check for assortativity 
    % the condition for PAM holds if σ_A < 1
    if param_vals(3) < 1
        assortativity = "possitive";
    else
        assortativity = "negative";
    end

    [dif_eqs, w, profits, vars, params, x_bounds, y_bounds, initial_conditions] = generate_model(assortativity);
    
	% Define an initial guess for firm size upper and lower bound
	firm_size_lower = 0;
	firm_size_upper = 1000000;
	initial_conditions(2) = (firm_size_lower + firm_size_upper)/2;
    if assortativity == "possitive"
        mu_last = y_bounds(2);
    else
        mu_last = y_bounds(1);
    end

    % Substitute parameters for their values in the differential equation
    for i = 1:numel(dif_eqs)
        dif_eqs(i) = subs(dif_eqs(i), params, param_vals);
    end

    % Solve the model for the first time (t = x, y = ( μ(x), θ(x) ) )
    [t,y] = solve_IVP(dif_eqs, vars, x_bounds, initial_conditions);

    n = 0 ;% Initialize iterator counter
    err = y(end, 1) - mu_last; % Initialize error
    
	while (abs(err) > tol) && (n < max_iter)
		
        n = n + 1;
		[t,y] = solve_IVP(dif_eqs, vars, x_bounds, initial_conditions);
		err =  y(end, 1) - mu_last;
		if err < 0
			firm_size_upper = initial_conditions(2);
			initial_conditions(2) = (firm_size_lower + firm_size_upper)/2;
		else
			firm_size_lower = initial_conditions(2);
			initial_conditions(2) = (firm_size_lower + firm_size_upper)/2;
		end
		fprintf( "Iteration: %d, err=%d, θ(1) = %d, μ(100) = %d \n",n, err, initial_conditions(2), y(end, 1) )
	end
    
    % Calculate profits and wages
    
    % Substitute parameters for their values in algebraic equation
    w = subs(w, params, param_vals);
    profits = subs(profits, params, param_vals);
    
    syms x % Re-generate x as symbolic
    wages = zeros(size(t));
    Pi = zeros(size(t));
    for i = 1:numel(t)
        wages(i) =  subs(w, [x vars], [t(i) y(i, :)]);
        Pi(i) =  subs(profits, [x vars], [t(i) y(i, :)]);
    end
    % Evaluate algebraic equations to obtain profits and wages
   
    y = [y, wages, Pi];

end