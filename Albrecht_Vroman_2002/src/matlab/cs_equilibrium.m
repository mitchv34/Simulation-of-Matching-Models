function cs_eq = cs_equilibrium(params)
%CS_EQUILIBRUM Summary of this function goes here
%   Detailed explanation goes here

    s_1   = params.s_1;
    s_2   = params.s_2;
    p     = params.p;
    b     = params.b;
    beta  = params.beta;
    delta = params.delta;
    c     = params.c;
    r     = params.r;

    
    syms theta gama; % naming gamma as gama to avoid confusion with matlab funciton `gamma`
    % Load arrival rate function 
    m = m_fun(theta);

    % Solve for the value of theta
    % define the equation that would give us theta

    eq_theta = c - ( m * (1 - beta )*(s_1 - c - b)) / ( theta * (r + delta + beta*m ) );
    %  solve for theta
    theta_eq = vpasolve(eq_theta == 0, theta, true);
    m_theta_eq = m(theta_eq);
    
    % solve for gamma
    %   Write the expression for phi in terms of gamma (using m_theta_eq)
    phi(gama) =( (1-gama)*p*m_theta_eq + (p - gama) * delta )/ (gama * m_theta_eq * (1 - p));
    
    % Write the equation that would give us gamma
    eq_gamma = gama*(r+delta)*(s_1-c-b)-(1-gama)*(s_2 - s_1)*(r+delta+m_theta_eq*phi(gama)*beta);
    
    gamma_eq = vpasolve(eq_gamma == 0, gama, true);
    % Eliminate the unadmisible solutions
    gamma_eq = gamma_eq( abs(gamma_eq) <= 1);
    phi_eq = phi(gamma_eq);

    % Solve for u in equilibrium
    u_eq = delta * (1 - p) / ((1 - gamma_eq)*(delta + m_theta_eq));
    
    % Value of unemployment for low skilled workers
    rUs_1 = (b*(r + delta) + m_theta_eq*beta*phi_eq*(s_1 - c) )/(r + delta + m_theta_eq*phi_eq*beta);
    % Value of unemployment for hihg-skilled workers
    rUs_2 = (b*(r+delta) + m_theta_eq * beta *(  phi_eq * s_1 + (1-phi_eq) * s_2 - c))/(r+delta+m_theta_eq*beta);
    
    % Calculate the value of wages 
    w_11 = wages(s_1, rUs_1, params);
    w_21 = wages(s_1, rUs_2, params);
    w_22 = wages(s_2, rUs_2, params);

    % TODO: Calculate Agregate Production

    % store equilibrium values in a struct
    cs_eq = struct('theta', theta_eq, 'm', m_theta_eq, 'u', u_eq, 'gamma', gamma_eq, 'phi', phi_eq, 'w_11', w_11, 'w_21', w_21, 'w_22', w_22);
end
