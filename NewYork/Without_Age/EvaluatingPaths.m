%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Evaluating Bootstrap predictions

TotalHospBoot = zeros(NumberOfAgeClasses,NSamples); 
TotalInfectionsBoot = zeros(NumberOfAgeClasses,NSamples);
TotalDeathsBoot = zeros(NumberOfAgeClasses,NSamples);

NewCasesBoot = zeros(length(t_actual),NSamples);
NewDeathsBoot = zeros(length(t_actual),NSamples);
NewHospBoot = zeros(length(t_actual),NSamples);

parfor ll = 1:NSamples
params2 = paramsOld;
p = params2.p;
Hosp = Hosp95(ll)*ones(length(t_actual),1);
% Hosp(1:end-4) = min(1,Hospitalization./trajectories(1:size(trajectories,1)-4,ll));
Hosp = min(20,Hosp/GetWorse_M);
params2.factorWorse = @(t)interp1(t_actual,Hosp,t);

% %%% Estimating the transmission constant parameters (M,H,I), the initial
% %%% infecve population (I_M0) and the transmission matrix:
yinitB = yinitOld;
I_M = PARAMS1Boot(ll,1);
yinitB(2*NumberOfAgeClasses+1:3*NumberOfAgeClasses) = (1-p)*I_M*PropInfections;
yinitB(3*NumberOfAgeClasses+1:4*NumberOfAgeClasses) = p*I_M*PropInfections;
yinitB(1:NumberOfAgeClasses) = Proportion-I_M*PropInfections;
params2.a = PARAMS1Boot(ll,2);
params2.b = PARAMS1Boot(ll,3);
params2.c = PARAMS1Boot(ll,4);

Beta = BETABoot(ll,:);
yinit2 = yinitB;
yb2 = zeros(length(t_actual),8*NumberOfAgeClasses);
yb2(1,:) = yinitB;
auxA = zeros(1,length(Beta)-1);
for jj =1:day-1
t_actual2 = t_actual(jj:jj+1);
beta2 = @(t)interp1(t_actual2,Beta(jj:jj+1),t);
[~,y2B] = ode45(@(t,y)seir_death_age_beta2(t,y,params2,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2B(end,:);
yb2(jj+1,:) = yinit2;
end

for ss = 1:nmax-1

if ss == nmax-1
aux = (day-1)+(ss-1)*dt+1:length(t_actual)-1;
aux2 = (day-1)+(ss-1)*dt+1:length(t_actual);    
else
aux = (day-1)+(ss-1)*dt+1:(day-1)+ss*dt;
aux2 = (day-1)+(ss-1)*dt+1:(day-1)+ss*dt+1;
end
params2.a = PARAMS2Boot(ss,1,ll);
params2.b = PARAMS2Boot(ss,2,ll);
params2.c = PARAMS2Boot(ss,3,ll);

for jj = aux
t_actual2 = t_actual(jj:jj+1);
beta2 = @(t)interp1(t_actual2,Beta(jj:jj+1),t);
[~,y2] = ode45(@(t,y)seir_death_age_beta2(t,y,params2,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2(end,:);
yb2(jj+1,:) = yinit2;
end
end
NewCasesBoot(:,ll) = ...
    sigma*sum(yb2(:,NumberOfAgeClasses+1:2*NumberOfAgeClasses),2);
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
      sum(sigma*yb(:,NumberOfAgeClasses+1:2*NumberOfAgeClasses))'*N;
TotalDeathsBoot(:,ll) = sum(Deathsb)'*N;
end