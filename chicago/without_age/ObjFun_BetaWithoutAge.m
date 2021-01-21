%Determine the beta rate taking into account the previous rates and without considering the age range.


function f = ObjFun_BetaWithoutAge(t_actual,params,data,options,priors,...
                                             yinit,unknowns)
Number = params.NumberOfAgeClasses;
params.a = unknowns(1);
p = params.p;

tspan = [t_actual(1),t_actual(end)];

N = params.N;


sigma = params.sigma;

[t,y]=ode45(@(t,y)seir_death_age_beta_b2(t,y, params),tspan,yinit,options);
NewInfections = p*sigma*sum(y(:,Number+1:2*Number),2)*N;
NewInfections = interp1(t,NewInfections,t_actual(2:end)');

% log-Poisson Misfit of Likelihood
Stirling = 0.5*log(2*pi*data(:,1)) + data(:,1).*log(data(:,1)) - data(:,1);
f = data(:,1).*log(NewInfections) - NewInfections - Stirling;
f = [f;1E-10*(unknowns-priors)'];
f(isnan(f)~=0)=zeros;