function  yprime = seir_Vaccination(t,y,params,beta)

%    dS/dt = -S*(beta_M*I_M + beta_H*I_H + beta_I*I_I)
%    dE/dt = S*(beta_M*I_M + beta_H*I_H + beta_I*I_I) - sigma*E 
%    dI_M/dt = sigma*E -(nu_M + mu_M + gamma_M)*I_M
%    dI_H/dt = gamma_M*I_M -(nu_H + mu_H + gamma_H)*I_H
%    dI_I/dt = gamma_H*I_H -(nu_I + mu_I)*I_I
%    dR/dt = nu_M*I_M + nu_H*I_H + nu_H*I_H
%    dD/dt = mu_M*I_M + mu_H*I_H + mu_H*I_H
%
%   yprime = [dS/dt  dE/dt dI_M/dt dI_H/dt dI_I/dt dR/dt dD/dt]'
%
%   input:
%       t  current time
%       y vector of current soln values
% y(1) = S, y(2) = E, y(3) = I_M, y(4) = I_H, y(5) = I_I, y(6) = R, 
% y(7) = D
%
%     parameters in "params" 
%    beta_M, beta_H, beta_I, N, sigma, gamma_M, gamma_H, nu_M, nu_H, nu_I, 
%    mu_M, mu_H, mu_I, R_zero_array (table of values) 
%
%   output: (col vector)
%    yprime(1) = dS/dt
%    yprime(2) = dE/dt
%    yprime(3) = dI_M/dt
%    yprime(4) = dI_H/dt
%    yprime(5) = dI_I/dt
%    yprime(6) = dR/dt
%    yprime(7) = dD/dt

%    beta_M, beta_H, beta_I, N, sigma, gamma_M, gamma_H, nu_M, nu_H, nu_I, 
%    mu_M, mu_H, mu_I,
%    R_zero_array (table of values)

%    R_zero_array = params.R_zero_array;
%     min_t = R_zero_array(1,1);
%       n_table = length( R_zero_array(:,1) );
%     max_t = R_zero_array(n_table,1);
%      t_val = max( min_t, min( t, max_t) );

%    R_zero = interp1( R_zero_array(:,1), R_zero_array(:,2), t_val);

a = params.a.*params.b;
b = params.a.*params.c;

factor = params.factorWorse;
factorD = params.factorDeath;

Number = params.NumberOfAgeClasses;
beta_M = beta(t)*params.beta_M;
beta_H = beta(t)*a*params.beta_H;
beta_I = beta(t)*b*params.beta_I;
% N = params.N/3;
sigma = params.sigma;
gamma_M = factor(t).*params.GetWorse_M;
gamma_H = params.GetWorse_H;
mu_M = params.Death_M;
mu_H = params.Death_H;
mu_I = factorD(t).*params.Death_I;
nu_M = ones-gamma_M;%params.Recovery_M;
nu_H = params.Recovery_H;
nu_A = params.Recovery_A;
nu_I = ones-mu_I;%params.Recovery_I;
p = params.p;
q = params.q;
beta_A = q*beta_M;
Vaccine = params.VaccinationRate;
nu = Vaccine(t);


%    beta_M = R_zero*gamma; 
S = y(1:Number);
% V = y(Number+1:2*Number);
E = y(2*Number+1:3*Number);
I_A = y(3*Number+1:4*Number);
I_M = y(4*Number+1:5*Number);
I_H = y(5*Number+1:6*Number);
I_I = y(6*Number+1:7*Number);
%      R = y(6);
%      D = y(7);
if sum(S) < 0.4
nu = zeros*nu;    
end
 yprime = zeros(9,1);
for jj = 1:Number
yprime(jj) = -S(jj)*(beta_M(jj,:)*I_M +beta_H(jj,:)*I_H +beta_I(jj,:)*I_I + beta_A(jj,:)*I_A) - nu(jj)*S(jj);
yprime(Number+jj) = nu(jj)*S(jj);
yprime(2*Number+jj) = S(jj)*(beta_M(jj,:)*I_M...
                  + beta_H(jj,:)*I_H + beta_I(jj,:)*I_I + beta_A(jj,:)*I_A) - sigma*E(jj);
yprime(3*Number+jj) = (1-p)*sigma*E(jj)-nu_A*I_A(jj);              
yprime(4*Number+jj) = p*sigma*E(jj)-(nu_M(jj)+mu_M(jj)+gamma_M(jj))*I_M(jj);
yprime(5*Number+jj) = ...
    gamma_M(jj)*I_M(jj) - (nu_H(jj)+mu_H(jj)+gamma_H(jj))*I_H(jj);
yprime(6*Number+jj) = ...
    gamma_H(jj)*I_H(jj)-(nu_I(jj) + mu_I(jj))*I_I(jj);
yprime(7*Number+jj) = nu_M(jj)*I_M(jj)+nu_H(jj)*I_H(jj)+nu_I(jj)*I_I(jj)+nu_A*I_A(jj);
yprime(8*Number+jj) = mu_M(jj)*I_M(jj)+mu_H(jj)*I_H(jj)+mu_I(jj)*I_I(jj);
end
