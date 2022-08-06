function es_eq = es_equilibrium(params)
    %ES_EQUILIBRUM Summary of this function goes here
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

    % Equation for phi in terms of theta and gamma
    phi(theta, gama) = (p * (1-gama) * m(theta) + delta*(p-gama)) / (m(theta)*(gama+p-2*gama*p));
    
    % Equation for u in terms of  theta and gamma
    u(theta,gama) = ( delta*(gama+p-2*gama*p)) / ( gama*(1-gama)*(m(theta)+ 2*delta));
    % System of equations to obtain equilibrium values of theta and gamma 
    eq_1 = (m(theta)/theta)*gama*( ((1-beta)*(s_1-c-b)) / (r+delta+m(theta)*phi(theta,gama)*beta)) - c;
    eq_2 = (m(theta)/theta)*(1-gama)*( ((1-beta)*(s_2-c-b)) / (r+delta+m(theta)*(1-phi(theta,gama))*beta)) - c;
    
    [theta_eq, gamma_eq] = vpasolve([eq_1 == 0, eq_2 == 0], [theta, gama], [0.3, 0.4]);

    m_theta_eq = m(theta_eq);
    phi_eq = phi(theta_eq, gamma_eq);
    u_eq = u(theta_eq, gamma_eq);


    % Value of unemployment for low skilled workers
    rUs_1 = (b*(r + delta) + m(theta_eq)*beta*phi_eq*(s_1 - c) )/(r + delta + m(theta_eq) * phi_eq * beta);
    % Value of unemployment for hihg-skilled workers
    rUs_2 = (b*(r + delta) + m(theta_eq)*beta*(1 - phi_eq)* (s_2 - c) )/( r + delta + m(theta_eq) * beta );
    
    % Calculate the value of wages 
    w_11 = wages(s_1, rUs_1, params);
    w_21 = NaN;
    w_22 = wages(s_2, rUs_2, params);

    % TODO: Calculate Agregate Production

    % store equilibrium values in a struct
    es_eq = struct( 'theta', theta_eq,...
                    'm', m_theta_eq,...
                    'u', u_eq,...ßßßßß
                    'gamma', gamma_eq,...
                    'phi', phi_eq,...
                    'w_11', w_11,...
                    'w_21', w_21,...
                    'w_22', w_22);
end % es_equilibrium