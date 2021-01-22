

function R0 = basic_reproduction_rate_beta(S,params,beta,t)

a = params.b;
b = params.c;

factor = params.factorWorse;
factorD = params.factorDeath;

Number = params.NumberOfAgeClasses;
beta_M = beta*params.beta_M;
beta_H = beta*a*params.beta_H;
beta_I = beta*b*params.beta_I;

sigma = params.sigma;
gamma_M = factor(t).*params.GetWorse_M;
gamma_H = params.GetWorse_H;
mu_M = params.Death_M;
mu_H = params.Death_H;
mu_I = factorD(t).*params.Death_I;
nu_M = ones-gamma_M;%params.Recovery_M;
nu_H = params.Recovery_H;
nu_I = ones-mu_I;%params.Recovery_I;
nu_A = params.Recovery_A;
q = params.q;
beta_A = q*beta_M;
p = params.p;

f = zeros(5*Number);
for jj = 1:Number
for ii =1:Number
f(jj,(ii-1)*5+1:ii*5) = [0,beta_A(jj,ii)*S(jj),beta_M(jj,ii)*S(jj),beta_H(jj,ii)*S(jj),beta_I(jj,ii)*S(jj)];
end
end
v = zeros(5*Number);

for jj = 1:Number
aux = diag([sigma,nu_A,nu_M(jj) + mu_M(jj) + gamma_M(jj),nu_H(jj) + mu_H(jj)...
                                       + gamma_H(jj),nu_I(jj) + mu_I(jj)]);
aux = aux - [zeros(1,5);diag([(1-p)*sigma,0,gamma_M(jj),gamma_H(jj)]),zeros(4,1)];
v((jj-1)*5+1:jj*5,(jj-1)*5+1:jj*5) = aux;
end
v(3,1) = -p*sigma;
aux = f/v;
eigen = eig(aux);
R0 = max(eigen);