function  yprime = seir_death_age_beta2(t,y,params,beta)


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

 
S = y(1:Number);
E = y(2*Number+1:3*Number);
I_A = y(3*Number+1:4*Number);
I_M = y(4*Number+1:5*Number);
I_H = y(5*Number+1:6*Number);
I_I = y(6*Number+1:7*Number);

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
