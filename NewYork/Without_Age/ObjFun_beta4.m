function f = ObjFun_beta4(t_actual,params,unknowns,data,options,priors,beta,yinit)
Number = params.NumberOfAgeClasses;
tspan = [t_actual(1),t_actual(end)];
beta = @(t)interp1(t_actual,[beta;unknowns],t);
p = params.p;
N = params.N;

sigma = params.sigma;

[~,y] = ode45(@(t,y)seir_death_age_beta2(t,y, params,beta),tspan,yinit,options);
NewInfections = p*sigma*sum(y(end,Number+1:2*Number),2)*N;

% f = (NewInfections-data(end,1));
% f = [f;10*(N*NewDeaths-data(end,2))];
% f = [f;1E-10*(unknowns2-priors)'];

f = data(end,1).*log(NewInfections) -NewInfections - sum(log(1:data(end,1)));
% f  = [f,10*(data(end,2).*log(NewDeaths) -NewDeaths - sum(log(1:data(end,2))))];
f = [f,1E-10*(unknowns-priors)];