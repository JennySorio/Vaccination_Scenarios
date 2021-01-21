function f = ObjFun_BetaM5(t_actual,params,data,options,priors,...
                                             yinit,unknowns,PropInfections)
Number = params.NumberOfAgeClasses;
params.a = unknowns(1);

bb = unknowns(2:end);
beta_M = 0.5*diag(PropInfections);
for jj = 1:Number-1
beta_M(jj,jj+1:end) = PropInfections(jj)*bb(1:end-jj+1);
end
beta_M = beta_M+beta_M';
params.beta_M = params.a.*beta_M;
params.beta_H = params.a.*beta_M;
params.beta_I = params.a.*beta_M;

tspan = [t_actual(1),t_actual(end)];

N = params.N;
p = params.p;

sigma = params.sigma;

[t,y]=ode45(@(t,y)seir_death_age_beta_b3(t,y, params),tspan,yinit,options);
NewInfections = p*sigma*sum(y(:,Number+1:2*Number),2)*N;
NewInfections = interp1(t,NewInfections,t_actual(2:end)');

% log-Poisson Misfit of Likelihood
Stirling = 0.5*log(2*pi*data(:,1)) + data(:,1).*log(data(:,1)) - data(:,1);
f = data(:,1).*log(NewInfections) - NewInfections - Stirling;
f = [f;1E-10*(unknowns-priors)'];
