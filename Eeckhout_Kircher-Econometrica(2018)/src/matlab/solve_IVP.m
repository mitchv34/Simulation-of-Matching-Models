% This function solves the IVP given an initial contions and parameter values
function [t,y] = solve_IVP(eqs, vars, x_bounds, initial_conditions)

    % Defining the ODEProblem object
    
    % Find the mass matrix M and vector of the right side F for this system
    [M, F] = massMatrixForm(eqs, vars);
    % To simplify further computations, rewrite the system in the form x'(t)=f(t,x(t)).
    f = M\F;
    % Convert f to a MATLAB function handle by using odeFunction.
    prob =  odeFunction(f,vars);
    % Solve using an ODE solver (in this case im using ode15s)
    [t,y] = ode15s(prob, x_bounds, initial_conditions);

end