%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bootstrap Sampling
NSamples = 200;

trajectories = zeros(length(t_actual)-1,NSamples);
for jj = 1:NSamples
trajectories(:,jj) = poissrnd(NewCases(2:end));
end

%% Estimating Bootstrap Parameters

PARAMS1Boot = zeros(NSamples,2);
BETABoot = zeros(NSamples,length(BETA));
R0Boot = zeros(NSamples,length(t_actual)-1);
AUX = zeros(size(BETABoot));
t_actualB = t_actual(day:end);

Death_M = params.Death_M;
Death_H = params.Death_H;
Death_I = params.Death_I;
GetWorse_M = params.GetWorse_M;
N = params.N;

TotalHospBoot = zeros(NumberOfAgeClasses,NSamples); 
TotalInfectionsBoot = zeros(NumberOfAgeClasses,NSamples);
TotalDeathsBoot = zeros(NumberOfAgeClasses,NSamples);

NewCasesBoot = zeros(length(t_actual),NSamples);
NewDeathsBoot = zeros(length(t_actual),NSamples);
NewHospBoot = zeros(length(t_actual),NSamples);

Hospitalization = data1(1:end-4,3);

dt = 20;
nmax = floor(length(t_actual)/dt);
PARAMS2Boot = zeros(nmax-1,NSamples);

parfor ll = 1:NSamples
params2=paramsOld;
p = params2.p;

%%%% Estimating the transmission constant parameters (M,H,I), the initial
%%%% infecve population (I_M0) and the transmission matrix:
yinitB = yinitOld;
LB = [1/N,1E-3];
UB = [50/N,10];
unknowns0=[1/N,beta];
priors = unknowns0;
t_actual2 = t_actual(1:day);

OF = @(unknowns)ObjFun_InitialPopBetaWithoutAge(t_actual2,params2,...
trajectories(1:day-1,ll),options,priors,yinitB,Proportion,PropInfections,unknowns);

unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);
I_M = unknowns(1);
yinitB(2*NumberOfAgeClasses+1:3*NumberOfAgeClasses) = (1-p)*I_M*PropInfections;
yinitB(3*NumberOfAgeClasses+1:4*NumberOfAgeClasses) = p*I_M*PropInfections;
yinitB(1:NumberOfAgeClasses) = Proportion-I_M*PropInfections;
params2.a = unknowns(2);

PARAMS1Boot(ll,:) = unknowns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Estimating a time-dependent transmission parameter (beta) from day 1
%%%% until 19:

Beta = zeros(length(t_actual),1);
Beta(1) = ones; 
unknowns0 = unknowns(2);
priors = unknowns0;
yinit2 = yinitB;
yb2 = zeros(length(t_actual),8*NumberOfAgeClasses);
yb2(1,:) = yinit;
auxA = zeros(1,length(Beta)-1);
for jj =1:day-1
t_actual2 = t_actual(jj:jj+1);
OF2 = @(unknowns2)ObjFun_beta4(t_actual2,params2,unknowns2,...
                    trajectories(jj,ll),options,priors,Beta(jj),yinit2);
unknowns = lsqnonlin(OF2,unknowns0,1E-2,10,options2);
unknowns0 = unknowns;
priors = unknowns;
Beta(jj+1) = unknowns;
beta2 = @(t)interp1(t_actual2,Beta(jj:jj+1),t);
[~,y2B] = ode45(@(t,y)seir_death_age_beta2(t,y,params2,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2B(end,:);
yb2(jj+1,:) = yinit2;
disp(num2str(unknowns))
auxA(jj) = ...
basic_reproduction_rate_beta(Proportion,params2,Beta(jj+1),t_actual(jj+1));
end


for ss = 1:nmax-1

LB = 1E-3;
UB = 10;
unknowns0= beta;
priors = unknowns0;

%%% Estimating the transmission constant parameters (M,H,I), the initial
%%% infecve population (I_M0) and the transmission matrix:
if ss == nmax-1
aux = (day-1)+(ss-1)*dt+1:length(t_actual)-1;
aux2 = (day-1)+(ss-1)*dt+1:length(t_actual);    
else
aux = (day-1)+(ss-1)*dt+1:(day-1)+ss*dt;
aux2 = (day-1)+(ss-1)*dt+1:(day-1)+ss*dt+1;
end
OF = @(unknowns)ObjFun_BetaWithoutAge(t_actual(aux2),params2,...
           trajectories(aux,ll),options,priors,yinit2,unknowns);
unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);
params2.a = unknowns(1);

PARAMS2Boot(ss,ll) = unknowns;

unknowns0 = ones;
priors = unknowns0;
for jj = aux
t_actual2 = t_actual(jj:jj+1);
OF2 = @(unknowns)ObjFun_beta4(t_actual2,params2,unknowns,...
                  trajectories(jj,ll),options,priors,Beta(jj),yinit2);
unknowns = lsqnonlin(OF2,unknowns0,1E-2,10,options3);
unknowns0 = unknowns;
priors = unknowns;
Beta(jj+1,:) = unknowns;
beta2 = @(t)interp1(t_actual2,Beta(jj:jj+1),t);
[~,y2] = ode45(@(t,y)seir_death_age_beta2(t,y,params2,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2(end,:);
yb2(jj+1,:) = yinit2;
disp(num2str(unknowns))
auxA(jj) = ...
basic_reproduction_rate_beta(Proportion,params2,Beta(jj+1),t_actual(jj+1));
end
end
AUX(ll,:) = Beta;
R0Boot(ll,:) = auxA;
NewCasesBoot(:,ll) = ...
    p*sigma*sum(yb2(:,NumberOfAgeClasses+1:2*NumberOfAgeClasses),2);

factor = zeros(length(t_actual),1);
factorD = factor;
factorW = params2.factorWorse;
factorDeath = params2.factorDeath;
for jj = 1:length(t_actual)
factor(jj) = factorW(t_actual(jj));
factorD(jj) = factorDeath(t_actual(jj));
end  

%%% Total Number of Deaths for each day
NewHospBb = zeros(size(yb(:,1),1),NumberOfAgeClasses);
Deathsb = zeros(length(yb(:,1)),NumberOfAgeClasses);

for jj=1:NumberOfAgeClasses
NewDeathsBoot(:,ll) = NewDeathsBoot(:,ll)...
+ Death_M(jj)*yb2(:,3*NumberOfAgeClasses+jj)...
+ Death_H(jj)*yb2(:,4*NumberOfAgeClasses+jj)...
+ Death_I(jj)*factorD.*yb2(:,5*NumberOfAgeClasses+jj);
NewHospBoot(:,ll) = NewHospBoot(:,ll)...
+ GetWorse_M(jj)*factor.*yb2(:,3*NumberOfAgeClasses+jj);

NewHospBb(:,jj) = GetWorse_M(jj)*factor.*yb2(:,3*NumberOfAgeClasses+jj);
Deathsb(:,jj) = Death_M(jj)*yb2(:,3*NumberOfAgeClasses+jj)...
+Death_H(jj)*yb2(:,4*NumberOfAgeClasses+jj)...
+Death_I(jj)*factorD.*yb2(:,5*NumberOfAgeClasses+jj);
end

%%% Total Hospitalizations for each Age Range
TotalHospBoot(:,ll) = sum(NewHospBb)'*N; 
TotalInfectionsBoot(:,ll) = ...
      sum(p*params.sigma*yb(:,NumberOfAgeClasses+1:2*NumberOfAgeClasses))'*N;
TotalDeathsBoot(:,ll) = sum(Deathsb)'*N;
end
% plot(t_actual(2:day),trajectories(1:day-1,1),'b',t_actual(2:day),NewCasesBoot(1:day-1,:)*N,'r')
BETABoot = AUX;
