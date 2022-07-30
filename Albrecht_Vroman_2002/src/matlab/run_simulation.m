% Define parameters
s_1 = 1;
s_2 = 1.30;
p = 2/3;
b = 0.1;
beta = 0.5;
delta = 0.2;
c = 0.3;
r = 0.05;
% Define m(theta)
% 

% Cross Skill Equilibrium

syms theta gamma;
syms m(theta) phi(gamma);

m(theta) = 2*theta^(1/2);

% Solve for the value of theta

% define the equation that would give us theta
eq_1 = c - ( m(theta)*(1 - beta)*(s_1 - c - b)) / ( theta * (r + delta + beta*m(theta)) ) ;

% solve for theta
theta_eq = vpasolve(eq_1 == 0, theta, true);
m_theta_eq = m(theta_eq);

% solve for gamma
% Write the expression for phi in terms of gamma (using m_theta_eq)
phi(gamma) =( (1-gamma)*p*m_theta_eq + (p - gamma) * delta )/ (gamma*m_theta_eq *(1-p));

% Write the equation that would give us gamma
eq_2 = gamma*(r+delta)*(s_1-c-b)-(1-gamma)*(s_2 - s_1)*(r+delta+m_theta_eq*phi(gamma)*beta);

gamma_eq = vpasolve(eq_2 == 0, gamma, [0,1]);
phi_eq = phi(gamma_eq);

% Solve for u in equilibrium
u_eq = delta * (1 - p) / ((1 - gamma_eq)*(delta + m_theta_eq));

%[theta_eq, m_theta_eq, u_eq, gamma_eq, phi_eq]

% Ex-post Segmentation Equilibrium

phi(gamma, theta) = (p*(1-gamma)*m(theta) + delta*(p-gamma))/( m(theta)*(gamma + p - 2*gamma*p) );

% Solve the system of equations that define gamma and theta
eq_3 =  (m(theta)/theta)*gamma*( ((1-beta)*(s_1-c-b)) / (r+delta+m(theta)*phi(theta,gamma)*beta)) - c;
eq_4 =  (m(theta)/theta)*(1-gamma)*( ((1-beta)*(s_2-c-b)) / (r+delta+m(theta)*(1-phi(theta,gamma))*beta)) - c;

sol = vpasolve([eq_3 == c, eq_4 == c], [theta, gamma]);
theta_eq = sol.theta;
gamma_eq = sol.gamma;
m_theta_eq = m(theta_eq);
phi_eq = phi(gamma_eq, theta_eq);


[theta_eq, gamma_eq]