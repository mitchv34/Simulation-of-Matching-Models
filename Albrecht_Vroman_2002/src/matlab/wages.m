function w = wages(y, rU, params)

    % Compute wages given skill required value of unemoployment and parameters
    w = params.beta*(y-params.c) + (1-params.beta) * rU;

end