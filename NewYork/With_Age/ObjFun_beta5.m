function f = ObjFun_beta5(t_actual,params,unknowns,data,options,priors,beta,yinit)
Number = params.NumberOfAgeClasses;
tspan = [t_actual(1),t_actual(end)];
beta = @(t)interp1(t_actual,[beta;unknowns],t);

N = params.N;

sigma = params.sigma;
p = params.p;

[~,y] = ode45(@(t,y)seir_death_age_beta3(t,y, params,beta),tspan,yinit,options);
NewInfections = p*sigma*sum(y(end,Number+1:2*Number),2)*N;


f = data(end,1).*log(NewInfections) -NewInfections - sum(log(1:data(end,1)));

f = [f,1E-1*(unknowns-priors)];
